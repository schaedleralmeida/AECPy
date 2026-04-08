"""
Módulo para definição do modelo da estrutura para análise
pelo Método da Rigidez Direta (MRD) no AECPy
"""

import numpy as np
import pandas as pd
from types import MappingProxyType
from typing import cast

from . import procedimentos as pmm
from .elemento import ElementoBase as Elemento, ElementoPP, ElementoPE, ElementoTP, ElementoTE, ElementoGR
from .no import NoBase, NoPP, NoPE, NoTP, NoTE, NoGR
from .secao import SecaoBase
from .unidades import Conversor
import matplotlib.pyplot as plt
from matplotlib.figure import Figure as MplFigure

from .graficos import (
    diagramas,
    modelo_2d,
    modelo_3d,
    tabela_rel,
    unidade_padrao_saida,
)


def numerar_gdls(nos: list[NoBase]) -> int:
    """
    Função para numerar os graus de liberdade dos nós
    Todos os graus de liberdade livres (F) são numerados de 1 a ngdlF
    Todos os graus de liberdade com deslocamentos definidos (S), sejam restritos (r) ou prescritos (p) são numerados de ngdlF+1 adiante
    """

    # variáveis auxiliares
    ii = 0
    all = list()

    # numerando os gdls livres
    for no in nos:
        igdl = list()
        gdlr = [no.idl(g) for g in no.deslocamentos_nulos]
        gdlp = [no.idl(g) for g in no.deslocamentos_prescritos]
        iS = sorted(gdlr + gdlp)
        for i in range(no.ngdl):
            if i not in iS:
                igdl.append(ii)
                ii += 1
            else:
                igdl.append(None)
        all.append(igdl)

    ngdlF = ii

    # numerando os gdls com deslocamentos definidos
    # e atribuindo a numeração à lista de nós
    for n in range(len(nos)):
        igdl = all[n]
        for i in range(len(igdl)):
            if igdl[i] is None:
                igdl[i] = ii
                ii += 1
        nos[n].igdl = igdl

    return ngdlF
    

def construir_SEL(
    nos: list[NoBase], els: list[Elemento], ngdl: int
) -> tuple[np.ndarray, np.ndarray]:
    """
    Função para construir a matriz de rigidez global e o vetor de forças nodais
    """

    K = np.zeros((ngdl, ngdl))
    F = np.zeros(ngdl)

    # Contribuição dos elementos à matriz de rigidez global e vetor de forças nodais
    for el in els:
        # matriz de rigidez do elemento em coordenadas locais
        K_el_local = el.Ke_local()

        # matriz de transformação de coordenadas global-local
        R = el.R()

        # cálculo da matriz de rigidez do elemento (no sistema global)
        # transformação da matriz de rigidez para os eixos globais
        K_el_global = pmm.transf_coord(K_el_local, R, inv=True)

        igdl = el.igdl
        # soma a contribuição do elememto à matriz de rigidez global
        pmm.espalhar(K_el_global, K , igdl)

        # Fornas nodais equivalentes em elementos carregados
        if el.carga or el.dTemp or el.inclui_peso_proprio:
            # Reações de engastamento perfeito (rep) em coordenadas locais
            rep_el_local = el.rep()
            # rep em coordenadas globais
            rep_el = pmm.transf_coord(rep_el_local, R, inv=True)
            # somando ao vetor de forças nodais equivalentes
            F[igdl] += -rep_el

    # Contribuição dos nós à matriz de rigidez global e vetor de forças nodais
    for no in nos:
        # nó com força externa
        if no.carga:
            F[no.igdl] += no.p
        # nó com deslocamento prescrito
        if no.deslocamentos_prescritos:
            for gdl, dp in no.deslocamentos_prescritos.items():
                ig = no.idg(gdl)
                F -= K[:, ig] * dp
        # nó com apoio elástico
        if no.apoio_elastico:
            for gdl, k in no.apoio_elastico.items():
                ig = no.idg(gdl)
                K[ig, ig] += k

    return K, F


def resolver_SEL(K: np.ndarray, F: np.ndarray, nos: list[NoBase], ngdlF: int) -> tuple[
    np.ndarray, np.ndarray
]:
    """
    Função para resolver o sistema de equações lineares
    """
    U = np.zeros(len(F))
    R = np.zeros(len(F))

    # Cálculo dos deslocamentos
    U[:ngdlF] = np.linalg.solve(K[:ngdlF, :ngdlF], F[:ngdlF])

    # Deslocamentos prescritos
    for no in nos:
        for gdl, dp in no.deslocamentos_prescritos.items():
            U[no.idg(gdl)] = dp

    # Reações
    R[ngdlF:] = K[ngdlF:, :ngdlF] @ U[:ngdlF] - F[ngdlF:]

    # Apoios elásticos
    for no in nos:
        for gdl, k_ae in no.apoio_elastico.items():
            ig = no.idg(gdl)
            R[ig] += -k_ae * U[ig]

    return U, R


def resultados_nos(nos: list[NoBase], U: np.ndarray, R: np.ndarray) -> list[pd.Series]:
    """
    Retorna lista de Series com deslocamentos e reações em cada nó.
    """
    deslocamentos = list(nos[0].gdls_globais)
    reacoes = ["R" + f for f in nos[0].forcas_globais]
    ngdl_no = nos[0].ngdl

    resultado = []
    for no in nos:
        igdl = no.igdl
        dados = np.zeros(2 * ngdl_no)
        dados[:ngdl_no] = U[igdl]
        dados[ngdl_no:] = R[igdl]
        resultado.append(pd.Series(dados, index=deslocamentos + reacoes))
    return resultado


def resultados_elementos(els: list[Elemento], U: np.ndarray, npts: int = 5) -> list[dict]:
    """
    Retorna lista de dicts com esforços internos em cada elemento.
    """
    resultados = []
    for el in els:
        x = np.linspace(0, el.L, npts)
        dl = el.dl(U)
        resultados.append(el.rel(dl, x))
    return resultados


# ---------------------------------------------------------------------------
# Funções auxiliares para leitura do arquivo Excel/ODS
# ---------------------------------------------------------------------------

def _parse_literal(val):
    """Converte o conteúdo de uma célula em um objeto Python.

    Exemplos de valores aceitos:
    - ``[0.0, 5.0]``          → lista de coordenadas
    - ``{'fx': 10.0}``        → dicionário de cargas
    - ``['ux', 'uz']``        → lista de graus de liberdade
    - Célula vazia / NaN     → ``None``
    """
    import ast
    try:
        if pd.isna(val):
            return None
    except (TypeError, ValueError):
        pass
    if isinstance(val, str):
        texto = val.strip()
        return ast.literal_eval(texto) if texto else None
    return val


def _parse_params(params_str) -> dict:
    """Converte uma string 'chave=valor, chave=valor' em dicionário.

    Exemplo: ``'A=0.02, I3=6.7e-5'`` → ``{'A': 0.02, 'I3': 6.7e-5}``
    Retorna um dicionário vazio quando o campo está em branco.
    """
    if params_str is None:
        return {}
    texto = str(params_str).strip()
    if not texto or texto.lower() in ("nan", "none"):
        return {}
    resultado: dict = {}
    for parte in texto.split(","):
        parte = parte.strip()
        if "=" in parte:
            chave, valor = parte.split("=", 1)
            resultado[chave.strip()] = float(valor.strip())
    return resultado


def _parse_bool(val) -> bool | None:
    """Converte o conteúdo de uma célula em booleano.

    Aceita ``True``/``False``, ``1``/``0``, ``'sim'``/``'yes'`` etc.
    Retorna ``None`` quando a célula está vazia.
    """
    try:
        if pd.isna(val):
            return None
    except (TypeError, ValueError):
        pass
    if isinstance(val, bool):
        return val
    if isinstance(val, (int, float)):
        return bool(int(val))
    return str(val).strip().lower() in ("true", "1", "sim", "yes")


class Modelo:
    """Modelo estrutural para análise pelo Método da Rigidez Direta (MRD).

    Parâmetros
    ----------
    tipo : str
        Tipo de modelo: 'PP' (pórtico plano), 'PE' (pórtico espacial),
        'TP' (treliça plana), 'TE' (treliça espacial) ou 'GR' (grelha).
    unidades : tuple[str, str], opcional
        Par (up_comprimento, up_forca) com as unidades de entrada dos dados do modelo.
        Exemplos: ('m', 'kN'), ('cm', 'kN'), ('mm', 'N').
        Quando informado, cria um conversor (`self.conv`) usado automaticamente
        nos métodos de visualização (tabela_nos, tabela_elemento, diagramas_elemento).
        Se omitido, nenhuma conversão de unidades é aplicada.
    """

    _tipos_no = {
        "PP": NoPP, "PE": NoPE, "TP": NoTP, "TE": NoTE, "GR": NoGR,
    }
    _tipos_elemento = {
        "PP": ElementoPP, "PE": ElementoPE, "TP": ElementoTP,
        "TE": ElementoTE, "GR": ElementoGR,
    }

    def __init__(self, tipo: str, unidades: tuple[str, str] | None = None) -> None:
        if tipo not in self._tipos_no:
            raise ValueError(
                f"Tipo de modelo inválido: '{tipo}'. "
                f"Use um dos: {list(self._tipos_no)}"
            )
        self.tipo = tipo
        if unidades is not None:
            self.conv: Conversor | None = Conversor(up_comprimento=unidades[0], up_forca=unidades[1])
        else:
            self.conv = None
        self._nos: dict[int, NoBase] = {}
        self._elementos: dict[int, Elemento] = {}
        self._secoes: dict[str, SecaoBase] = {}
        self._nos_lista: list[NoBase] | None = None
        self._elementos_lista: list[Elemento] | None = None
        # --- síntese (preenchidas em concluir) ---
        self.concluido: bool = False
        self.ngdlF: int | None = None
        # --- resultados (preenchidos em analisar) ---
        self.K: np.ndarray | None = None
        self.F: np.ndarray | None = None
        self.U: np.ndarray | None = None
        self.R: np.ndarray | None = None
        self._res_nos: dict[int, pd.Series] | None = None
        self._res_elementos: dict[int, dict] | None = None

    # ------------------------------------------------------------------
    # Adição de componentes

    def adicionar_no(
        self,
        id: int,
        coor,
        *,
        carga=None,
        deslocamentos_nulos=None,
        deslocamentos_prescritos=None,
        apoio_elastico=None,
    ) -> NoBase:
        """Cria e registra um nó no modelo.

        Parâmetros
        ----------
        id : int
            Identificador externo do nó.
        coor : array-like
            Coordenadas do nó.
        carga : dict ou array-like, opcional
            Cargas nodais (ex.: {'fx': 10.0, 'fz': -5.0}).
        deslocamentos_nulos : list, opcional
            GDLs com apoio (ex.: ['ux', 'uz']).
        deslocamentos_prescritos : dict, opcional
            GDLs com deslocamento prescrito (ex.: {'uz': 0.002}).
        apoio_elastico : dict, opcional
            GDLs com apoio elástico e rigidez (ex.: {'uz': 1e5}).
        """
        cls_no = self._tipos_no[self.tipo]
        no = cls_no(coor)
        if carga is not None:
            no.carga = carga
        if deslocamentos_nulos is not None:
            no.deslocamentos_nulos = deslocamentos_nulos
        if deslocamentos_prescritos is not None:
            no.deslocamentos_prescritos = deslocamentos_prescritos
        if apoio_elastico is not None:
            no.apoio_elastico = apoio_elastico
        self._nos[id] = no
        self.abrir()
        return no

    def adicionar_elemento(
        self,
        id: int,
        id_noI: int,
        id_noJ: int,
        sec,
        *,
        carga=None,
        dTemp=None,
        inclui_peso_proprio=None,
    ) -> Elemento:
        """Cria e registra um elemento no modelo.

        Parâmetros
        ----------
        id : int
            Identificador externo do elemento.
        id_noI : int
            ID do nó inicial (já adicionado ao modelo).
        id_noJ : int
            ID do nó final (já adicionado ao modelo).
        sec : Secao
            Seção transversal do elemento.
        carga : dict, opcional
            Cargas distribuídas (ex.: {'w2': -5.0}).
        dTemp : dict, opcional
            Variações de temperatura (ex.: {'dT0': 20.0}).
        inclui_peso_proprio : bool, opcional
            Se True, o peso próprio é incluído no carregamento.
        """
        cls_el = self._tipos_elemento[self.tipo]
        el = cls_el(self._nos[id_noI], self._nos[id_noJ], sec)
        self._registrar_secao(sec)
        if carga is not None:
            el.carga = carga
        if dTemp is not None:
            el.dTemp = dTemp
        if inclui_peso_proprio is not None:
            el.inclui_peso_proprio = inclui_peso_proprio
        self._elementos[id] = el
        self.abrir()
        return el

    def _registrar_secao(self, sec: SecaoBase) -> None:
        """Registra a seção no dicionário interno, evitando duplicatas."""
        if any(s is sec for s in self._secoes.values()):
            return
        if sec.nome != "":
            self._secoes[sec.nome] = sec
        else:
            i = 1
            while f"sec_{i}" in self._secoes:
                i += 1
            self._secoes[f"sec_{i}"] = sec

    # ------------------------------------------------------------------
    # Carregamento de arquivo Excel / ODS

    def carregar_excel(self, caminho: str) -> None:
        """Carrega os dados do modelo a partir de um arquivo Excel ou ODS.

        O arquivo deve conter as seguintes abas:

        **materiais** — colunas: ``label``, ``E`` e, opcionalmente, ``G``,
        ``pe``, ``cdt``, ``nome``.

        **secoes** — colunas: ``label``, ``tipo`` (nome da classe: ``SecaoBase``,
        ``Secao`` ou ``SecaoRetangular``), ``material`` (label do material,
        obrigatório para ``Secao`` e ``SecaoRetangular``), ``params`` (string
        com argumentos adicionais no formato ``chave=valor, chave=valor``,
        p. ex. ``A=0.02, I3=6.7e-5``).

        **nos** — colunas: ``id``, ``coor`` e, opcionalmente, ``carga``,
        ``deslocamentos_nulos``, ``deslocamentos_prescritos``, ``apoio_elastico``.
        Os valores de ``coor`` e das colunas opcionais devem ser escritos como
        literais Python (lista ou dicionário), p. ex. ``[0.0, 5.0]`` ou
        ``{'fx': 10.0, 'fz': -5.0}``.

        **elementos** — colunas: ``id``, ``id_noI``, ``id_noJ``, ``sec``
        (label da seção) e, opcionalmente, ``carga``, ``dTemp``,
        ``inclui_peso_proprio``.

        Parâmetros
        ----------
        caminho : str
            Caminho para o arquivo (``.xlsx``, ``.ods``, ``.xls``, …).
            Requer ``openpyxl`` para ``.xlsx`` e ``odfpy`` para ``.ods``.
        """
        xls = pd.ExcelFile(caminho)
        materiais = self._ler_materiais(xls)
        secoes = self._ler_secoes(xls, materiais)
        self._ler_nos(xls)
        self._ler_elementos(xls, secoes)

    def _ler_materiais(self, xls: pd.ExcelFile) -> dict:
        """Lê a aba 'materiais' e retorna um dicionário {label: Material}.

        Colunas esperadas: ``label``, ``E`` (obrigatórias).
        Colunas opcionais: ``G``, ``pe``, ``cdt``, ``nome``.
        """
        from .secao import Material

        if "materiais" not in xls.sheet_names:
            raise ValueError(
                "Aba 'materiais' não encontrada no arquivo. "
                f"Abas disponíveis: {xls.sheet_names}"
            )
        materiais: dict = {}

        df = pd.read_excel(xls, sheet_name="materiais", dtype=str)
        df = df.where(df.notna(), other=None)

        for _, linha in df.iterrows():
            label = str(linha["label"])

            E = float(linha["E"])
            G   = float(linha["G"])   if "G"   in linha.index and linha["G"]   is not None else 0.0
            pe  = float(linha["pe"])  if "pe"  in linha.index and linha["pe"]  is not None else 0.0
            cdt = float(linha["cdt"]) if "cdt" in linha.index and linha["cdt"] is not None else 0.0
            nome = str(linha["nome"]) if "nome" in linha.index and linha["nome"] is not None else ""

            materiais[label] = Material(E=E, G=G, pe=pe, cdt=cdt, nome=nome)

        return materiais

    def _ler_secoes(self, xls: pd.ExcelFile, materiais: dict) -> dict:
        """Lê a aba 'secoes' e retorna um dicionário {label: SecaoBase}.

        Colunas esperadas: ``label``, ``tipo``, ``params``.
        Coluna ``material`` é obrigatória para os tipos ``Secao`` e
        ``SecaoRetangular``.

        O campo ``tipo`` deve conter o nome exato da classe:
        ``SecaoBase``, ``Secao`` ou ``SecaoRetangular``.

        O campo ``params`` deve conter os demais argumentos no formato
        ``chave=valor, chave=valor`` (ex.: ``A=0.02, I3=6.7e-5``).
        """
        from .secao import Secao, SecaoRetangular

        _tipos_secao: dict[str, type] = {
            "SecaoBase": SecaoBase,
            "Secao": Secao,
            "SecaoRetangular": SecaoRetangular,
        }

        if "secoes" not in xls.sheet_names:
            raise ValueError(
                "Aba 'secoes' não encontrada no arquivo. "
                f"Abas disponíveis: {xls.sheet_names}"
            )
        secoes: dict = {}

        df = pd.read_excel(xls, sheet_name="secoes", dtype=str)
        df = df.where(df.notna(), other=None)

        for _, linha in df.iterrows():
            label = str(linha["label"])
            tipo_nome = str(linha["tipo"]).strip()

            cls = _tipos_secao.get(tipo_nome)
            if cls is None:
                raise ValueError(
                    f"Tipo de seção desconhecido: '{tipo_nome}'. "
                    f"Use um de: {list(_tipos_secao)}"
                )

            # Converte 'A=0.02, I3=6.7e-5' em {'A': 0.02, 'I3': 6.7e-5}
            params_str = linha["params"] if "params" in linha.index else None
            kwargs = _parse_params(params_str)
            kwargs["nome"] = label

            # Tipos que precisam de um objeto Material
            if issubclass(cls, Secao):
                mat_label = (
                    str(linha["material"]).strip()
                    if "material" in linha.index and linha["material"] is not None
                    else None
                )
                if mat_label is None:
                    raise ValueError(
                        f"A seção '{label}' (tipo '{tipo_nome}') requer a coluna 'material'."
                    )
                if mat_label not in materiais:
                    raise KeyError(
                        f"Material '{mat_label}' não encontrado para a seção '{label}'. "
                        f"Materiais disponíveis: {list(materiais)}"
                    )
                secoes[label] = cls(materiais[mat_label], **kwargs)
            else:
                secoes[label] = cls(**kwargs)

        return secoes

    def _ler_nos(self, xls: pd.ExcelFile) -> None:
        """Lê a aba 'nos' e adiciona cada nó ao modelo.

        Colunas obrigatórias: ``id`` e as colunas de coordenadas com os nomes
        dos eixos globais do tipo de nó (p. ex. ``x`` e ``z`` para pórtico
        plano; ``x``, ``y`` e ``z`` para pórtico espacial).

        Colunas opcionais: ``carga``, ``deslocamentos_nulos``,
        ``deslocamentos_prescritos``, ``apoio_elastico``.
        Esses campos devem ser escritos como literais Python,
        por exemplo ``{'fx': 10.0}`` ou ``['ux', 'uz']``.
        """
        if "nos" not in xls.sheet_names:
            raise ValueError(
                "Aba 'nos' não encontrada no arquivo. "
                f"Abas disponíveis: {xls.sheet_names}"
            )

        df = pd.read_excel(xls, sheet_name="nos")

        # Eixos globais do tipo de nó deste modelo (ex.: ('x', 'z') para PP)
        eixos = self._tipos_no[self.tipo].eixos_globais

        # Verifica se todas as colunas de coordenadas estão presentes
        for eixo in eixos:
            if eixo not in df.columns:
                raise ValueError(
                    f"Coluna '{eixo}' não encontrada na aba 'nos'. "
                    f"Para o tipo '{self.tipo}', as colunas de coordenadas esperadas são: {list(eixos)}"
                )

        for _, linha in df.iterrows():
            id_no = int(linha["id"])

            # Monta a lista de coordenadas na ordem correta dos eixos globais
            coor = [float(linha[eixo]) for eixo in eixos]

            # Campos opcionais: só são incluídos se a coluna existir e não estiver vazia
            opcoes = {}
            if "carga"                    in linha.index: opcoes["carga"]                    = _parse_literal(linha["carga"])
            if "deslocamentos_nulos"      in linha.index: opcoes["deslocamentos_nulos"]      = _parse_literal(linha["deslocamentos_nulos"])
            if "deslocamentos_prescritos" in linha.index: opcoes["deslocamentos_prescritos"] = _parse_literal(linha["deslocamentos_prescritos"])
            if "apoio_elastico"           in linha.index: opcoes["apoio_elastico"]           = _parse_literal(linha["apoio_elastico"])
            # Remove entradas cuja célula estava vazia (None)
            opcoes = {chave: valor for chave, valor in opcoes.items() if valor is not None}

            self.adicionar_no(id_no, coor, **opcoes)

    def _ler_elementos(self, xls: pd.ExcelFile, secoes: dict) -> None:
        """Lê a aba 'elementos' e adiciona cada elemento ao modelo.

        Colunas esperadas: ``id``, ``id_noI``, ``id_noJ``, ``sec``
        (label da seção, conforme definido na aba 'secoes').
        Colunas opcionais: ``carga``, ``dTemp``, ``inclui_peso_proprio``.
        """
        if "elementos" not in xls.sheet_names:
            raise ValueError(
                "Aba 'elementos' não encontrada no arquivo. "
                f"Abas disponíveis: {xls.sheet_names}"
            )

        df = pd.read_excel(xls, sheet_name="elementos")

        for _, linha in df.iterrows():
            id_el = int(linha["id"])
            id_noI = int(linha["id_noI"])
            id_noJ = int(linha["id_noJ"])
            sec_label = str(linha["sec"]).strip()

            if sec_label not in secoes:
                raise KeyError(
                    f"Seção '{sec_label}' não encontrada para o elemento {id_el}. "
                    f"Seções disponíveis: {list(secoes)}"
                )

            kwargs: dict = {}
            for coluna in ("carga", "dTemp"):
                if coluna in linha.index:
                    valor = _parse_literal(linha[coluna])
                    if valor is not None:
                        kwargs[coluna] = valor

            if "inclui_peso_proprio" in linha.index:
                kwargs["inclui_peso_proprio"] = _parse_bool(linha["inclui_peso_proprio"])

            self.adicionar_elemento(id_el, id_noI, id_noJ, secoes[sec_label], **kwargs)

    def abrir(self) -> None:
        """Reabre o modelo para edição, zerando todos os dados derivados.

        Chamado automaticamente ao adicionar nós ou elementos.
        """
        self.concluido = False
        self.ngdlF = None
        self._nos_lista = None
        self._elementos_lista = None
        self.K = self.F = self.U = self.R = None
        self._res_nos = self._res_elementos = None

    def concluir(self) -> None:
        """Conclui a montagem do modelo.

        Numera os graus de liberdade dos nós, congela as listas internas
        e calcula as variáveis de síntese do modelo.
        Deve ser chamado após adicionar todos os nós e elementos.
        """
        self._nos_lista = list(self._nos.values())
        self._elementos_lista = list(self._elementos.values())
        self.ngdlF = numerar_gdls(self._nos_lista)
        self.K = self.F = self.U = self.R = None
        self._res_nos = self._res_elementos = None
        self.concluido = True
        assert self._nos_lista is not None
        assert self._elementos_lista is not None
        assert self.ngdlF is not None

    def construir_SEL(self) -> None:
        """Constrói e armazena a matriz de rigidez global (K) e o vetor de forças (F).

        Requer que concluir() já tenha sido chamado.
        Os resultados ficam disponíveis em self.K e self.F.
        """
        if not self.concluido:
            self.concluir()
        self.K, self.F = construir_SEL(
            self._nos_lista, self._elementos_lista, self.ngdl  # type: ignore[arg-type]
        )

    def analisar(self, npts: int = 5) -> None:
        """Executa a análise completa da estrutura.

        Se o modelo ainda não foi concluído, chama concluir() automaticamente.
        Os resultados ficam disponíveis em res_nos e res_elementos.

        Parâmetros
        ----------
        npts : int, opcional
            Número de pontos por elemento para o cálculo dos esforços internos (default 5).
        """
        self.construir_SEL()
        self.U, self.R = resolver_SEL(
            self.K, self.F, self._nos_lista, cast(int, self.ngdlF)  # type: ignore[arg-type]
        )
        res_nos_list = resultados_nos(self._nos_lista, self.U, self.R)  # type: ignore[arg-type]
        res_els_list = resultados_elementos(self._elementos_lista, self.U, npts)  # type: ignore[arg-type]
        self._res_nos = dict(zip(self._nos.keys(), res_nos_list))
        self._res_elementos = dict(zip(self._elementos.keys(), res_els_list))

    # ------------------------------------------------------------------
    # Acesso

    @property
    def res_nos(self) -> dict[int, pd.Series] | None:
        """Resultados nos nós. None se o modelo não foi analisado."""
        return self._res_nos

    @property
    def res_elementos(self) -> dict[int, dict] | None:
        """Resultados nos elementos. None se o modelo não foi analisado."""
        return self._res_elementos

    @property
    def no(self) -> MappingProxyType:
        """Dicionário somente-leitura de nós. Acesso: modelo.no[id]"""
        return MappingProxyType(self._nos)

    @property
    def elemento(self) -> MappingProxyType:
        """Dicionário somente-leitura de elementos. Acesso: modelo.elemento[id]"""
        return MappingProxyType(self._elementos)

    @property
    def secoes(self) -> MappingProxyType:
        """Dicionário somente-leitura de seções registradas, indexado pelo nome."""
        return MappingProxyType(self._secoes)

    @property
    def n_nos(self) -> int:
        """Número de nós do modelo."""
        return len(self._nos)

    @property
    def n_elementos(self) -> int:
        """Número de elementos do modelo."""
        return len(self._elementos)

    @property
    def ngdl(self) -> int:
        """Total de graus de liberdade do modelo."""
        return len(self._nos) * self._tipos_no[self.tipo].ngdl

    @property
    def ngdlS(self) -> int | None:
        """Graus de liberdade impedidos (restritos + prescritos). None antes de concluir()."""
        if self.ngdlF is None:
            return None
        return self.ngdl - self.ngdlF

    # ------------------------------------------------------------------
    # Visualização

    def visualizar(self) -> MplFigure:
        """Cria visualização geométrica 2D ou 3D do modelo.

        Funciona em qualquer estado (não requer análise).
        """
        if not self._nos:
            raise RuntimeError("O modelo não possui nós.")
        nos = list(self._nos.values())
        els = list(self._elementos.values())
        if self._tipos_no[self.tipo].ndim == 3:
            return modelo_3d(nos, els)
        return modelo_2d(nos, els)

    def diagramas_elemento(
        self,
        id: int,
        resultados: list[str] | None = None,
        conv: Conversor | None = None,
    ) -> MplFigure:
        """Cria diagramas de esforços e deslocamentos ao longo de um elemento.

        Parâmetros
        ----------
        id : int
            Identificador externo do elemento.
        resultados : list[str], opcional
            Grandezas a plotar (ex.: ['N', 'V2', 'M3']). Se omitido, plota todas.
        conv : Conversor, opcional
            Conversor de unidades para exibição. Se omitido, usa o conversor do modelo.
        """
        if self._res_elementos is None:
            raise RuntimeError("Modelo não foi analisado. Chame analisar() primeiro.")
        if id not in self._res_elementos:
            raise KeyError(f"Elemento {id} não encontrado.")
        return diagramas(self._res_elementos[id], resultados, conv if conv is not None else self.conv, eltag=id)

    def tabela_nos(
        self,
        ids: list[int] | None = None,
        conv: Conversor | None = None,
    ) -> pd.DataFrame:
        """Retorna DataFrame com deslocamentos e reações nos nós.

        O índice do DataFrame são os IDs externos dos nós.

        Parâmetros
        ----------
        ids : list[int], opcional
            IDs dos nós a incluir. Se omitido, inclui todos os nós.
        conv : Conversor, opcional
            Conversor de unidades para exibição. Se omitido, usa o conversor do modelo.
        """
        if self._res_nos is None:
            raise RuntimeError("Modelo não foi analisado. Chame analisar() primeiro.")
        if ids is not None:
            dados = {id: self._res_nos[id] for id in ids}
        else:
            dados = self._res_nos
        df = pd.DataFrame(dados).T
        df.index.name = "nó"
        conv = conv if conv is not None else self.conv
        if conv is not None:
            rename: dict[str, str] = {}
            for col in df.columns:
                lookup = col[1:] if col.startswith("R") else col
                ups = unidade_padrao_saida(lookup)
                if ups:
                    try:
                        c = conv.para(ups)
                        df[col] = df[col] * c
                        rename[col] = f"{col} ({ups})"
                    except Exception:
                        pass
            df = df.rename(columns=rename)
        return df

    def tabela_elemento(
        self,
        id: int,
        resultados: list[str] | None = None,
        conv: Conversor | None = None,
    ) -> pd.DataFrame:
        """Retorna DataFrame com resultados ao longo de um elemento.

        Parâmetros
        ----------
        id : int
            Identificador externo do elemento.
        resultados : list[str], opcional
            Grandezas a incluir (ex.: ['x', 'N', 'M3']). Se omitido, inclui todas.
        conv : Conversor, opcional
            Conversor de unidades para exibição. Se omitido, usa o conversor do modelo.
        """
        if self._res_elementos is None:
            raise RuntimeError("Modelo não foi analisado. Chame analisar() primeiro.")
        if id not in self._res_elementos:
            raise KeyError(f"Elemento {id} não encontrado.")
        return tabela_rel(self._res_elementos[id], resultados or [], conv if conv is not None else self.conv)

    def __str__(self) -> str:
        linhas = [f"Modelo  tipo: {self.tipo}"]
        linhas.append(f"  status    : {'concluído' if self.concluido else 'aberto'}")
        linhas.append(f"  nós       : {len(self._nos)}")
        linhas.append(f"  elementos : {len(self._elementos)}")
        linhas.append(f"  seções    : {list(self._secoes)}")
        if self.concluido:
            linhas.append(f"  GDLs      : {self.ngdl}  (livres: {self.ngdlF}  impedidos: {self.ngdlS})")
            analisado = self._res_nos is not None
            linhas.append(f"  resultados: {'disponíveis' if analisado else 'não analisado'}")
        return "\n".join(linhas)
