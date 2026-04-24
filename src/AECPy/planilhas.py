"""
Módulo para leitura de dados estruturais a partir de planilhas Excel/ODS.

As funções ``ler_*`` recebem um :class:`pandas.DataFrame` já carregado e o
objeto :class:`~AECPy.modelo.Modelo` a ser populado, tornando-as independentes
do formato do arquivo de origem.  A função :func:`carregar_excel` é uma camada
de conveniência que abre o arquivo, lê cada aba e delega para as funções
específicas.
"""

from __future__ import annotations

import ast
from typing import TYPE_CHECKING

import pandas as pd

from .secao import SecaoBase

if TYPE_CHECKING:
    from .modelo import Modelo


# ---------------------------------------------------------------------------
# Helpers de conversão de célula
# ---------------------------------------------------------------------------

def _parse_literal(val):
    """Converte o conteúdo de uma célula em um objeto Python.

    Exemplos de valores aceitos:

    - ``[0.0, 5.0]``       → lista de coordenadas
    - ``{'fx': 10.0}``     → dicionário de cargas
    - ``['ux', 'uz']``     → lista de graus de liberdade
    - Célula vazia / NaN  → ``None``
    """
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
    """Converte uma string ``'chave=valor, chave=valor'`` em dicionário.

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


# ---------------------------------------------------------------------------
# Leitores de cada aba
# ---------------------------------------------------------------------------

def ler_materiais(df: pd.DataFrame) -> dict:
    """Constrói um dicionário ``{label: Material}`` a partir de um DataFrame.

    Colunas obrigatórias: ``label``, ``E``.
    Colunas opcionais: ``G``, ``pe``, ``cdt``, ``nome``.

    Parâmetros
    ----------
    df : pd.DataFrame
        DataFrame com os dados dos materiais (tipicamente a aba *materiais*
        do arquivo Excel/ODS), lido com ``dtype=str``.

    Retorna
    -------
    dict
        Dicionário ``{label: Material}``.
    """
    from .secao import Material

    df = df.where(df.notna(), other=None)
    materiais: dict = {}

    for _, linha in df.iterrows():
        label = str(linha["label"])
        E = float(linha["E"])
        G   = float(linha["G"])   if "G"   in linha.index and linha["G"]   is not None else 0.0
        pe  = float(linha["pe"])  if "pe"  in linha.index and linha["pe"]  is not None else 0.0
        cdt = float(linha["cdt"]) if "cdt" in linha.index and linha["cdt"] is not None else 0.0
        nome = str(linha["nome"]) if "nome" in linha.index and linha["nome"] is not None else ""
        materiais[label] = Material(E=E, G=G, pe=pe, cdt=cdt, nome=nome)

    return materiais


def ler_secoes(df: pd.DataFrame, materiais: dict) -> dict:
    """Constrói um dicionário ``{label: SecaoBase}`` a partir de um DataFrame.

    Colunas obrigatórias: ``label``, ``tipo``, ``params``.
    Coluna ``material`` é obrigatória para os tipos ``Secao`` e
    ``SecaoRetangular``.

    O campo ``tipo`` deve conter o nome exato da classe:
    ``SecaoBase``, ``Secao`` ou ``SecaoRetangular``.

    O campo ``params`` deve conter os argumentos adicionais no formato
    ``chave=valor, chave=valor`` (ex.: ``A=0.02, I3=6.7e-5``).

    Parâmetros
    ----------
    df : pd.DataFrame
        DataFrame com os dados das seções (tipicamente a aba *secoes*),
        lido com ``dtype=str``.
    materiais : dict
        Dicionário ``{label: Material}`` retornado por :func:`ler_materiais`.

    Retorna
    -------
    dict
        Dicionário ``{label: SecaoBase}``.
    """
    from .secao import Secao, SecaoRetangular

    _tipos_secao: dict[str, type] = {
        "SecaoBase": SecaoBase,
        "Secao": Secao,
        "SecaoRetangular": SecaoRetangular,
    }

    df = df.where(df.notna(), other=None)
    secoes: dict = {}

    for _, linha in df.iterrows():
        label = str(linha["label"])
        tipo_nome = str(linha["tipo"]).strip()

        cls = _tipos_secao.get(tipo_nome)
        if cls is None:
            raise ValueError(
                f"Tipo de seção desconhecido: '{tipo_nome}'. "
                f"Use um de: {list(_tipos_secao)}"
            )

        params_str = linha["params"] if "params" in linha.index else None
        kwargs = _parse_params(params_str)
        kwargs["nome"] = label

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


def ler_nos(df: pd.DataFrame, modelo: Modelo) -> None:
    """Lê nós a partir de um DataFrame e os adiciona ao modelo.

    Colunas obrigatórias: ``id`` e as colunas de coordenadas com os nomes dos
    eixos globais do tipo de nó (ex.: ``x`` e ``z`` para pórtico plano;
    ``x``, ``y`` e ``z`` para pórtico espacial).

    Colunas opcionais: ``carga``, ``deslocamentos_nulos``,
    ``deslocamentos_prescritos``, ``apoio_elastico``.
    Esses campos devem ser escritos como literais Python,
    por exemplo ``{'fx': 10.0}`` ou ``['ux', 'uz']``.

    Parâmetros
    ----------
    df : pd.DataFrame
        DataFrame com os dados dos nós (tipicamente a aba *nos*).
    modelo : Modelo
        Instância do modelo que receberá os nós.
    """
    eixos = modelo._tipos_no[modelo.tipo].eixos_globais

    for eixo in eixos:
        if eixo not in df.columns:
            raise ValueError(
                f"Coluna '{eixo}' não encontrada no DataFrame de nós. "
                f"Para o tipo '{modelo.tipo}', as colunas de coordenadas esperadas são: {list(eixos)}"
            )

    for _, linha in df.iterrows():
        id_no = int(linha["id"])
        coor = [float(linha[eixo]) for eixo in eixos]

        opcoes = {}
        if "carga"                    in linha.index: opcoes["carga"]                    = _parse_literal(linha["carga"])
        if "deslocamentos_nulos"      in linha.index: opcoes["deslocamentos_nulos"]      = _parse_literal(linha["deslocamentos_nulos"])
        if "deslocamentos_prescritos" in linha.index: opcoes["deslocamentos_prescritos"] = _parse_literal(linha["deslocamentos_prescritos"])
        if "apoio_elastico"           in linha.index: opcoes["apoio_elastico"]           = _parse_literal(linha["apoio_elastico"])
        opcoes = {chave: valor for chave, valor in opcoes.items() if valor is not None}

        modelo.adicionar_no(id_no, coor, **opcoes)


def ler_elementos(df: pd.DataFrame, secoes: dict, modelo: Modelo) -> None:
    """Lê elementos a partir de um DataFrame e os adiciona ao modelo.

    Colunas obrigatórias: ``id``, ``id_noI``, ``id_noJ``, ``sec``
    (label da seção, conforme definido em ``secoes``).
    Colunas opcionais: ``carga``, ``dTemp``, ``inclui_peso_proprio``.

    Parâmetros
    ----------
    df : pd.DataFrame
        DataFrame com os dados dos elementos (tipicamente a aba *elementos*).
    secoes : dict
        Dicionário ``{label: SecaoBase}`` retornado por :func:`ler_secoes`.
    modelo : Modelo
        Instância do modelo que receberá os elementos.
    """
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

        modelo.adicionar_elemento(id_el, id_noI, id_noJ, secoes[sec_label], **kwargs)


# ---------------------------------------------------------------------------
# Função de conveniência
# ---------------------------------------------------------------------------

def carregar_excel(caminho: str, modelo: Modelo) -> None:
    """Carrega os dados do modelo a partir de um arquivo Excel ou ODS.

    Lê as abas *materiais*, *secoes*, *nos* e *elementos* e popula o modelo
    chamando :func:`ler_materiais`, :func:`ler_secoes`, :func:`ler_nos` e
    :func:`ler_elementos`.

    O arquivo deve conter as seguintes abas:

    **materiais** — colunas: ``label``, ``E`` e, opcionalmente, ``G``,
    ``pe``, ``cdt``, ``nome``.

    **secoes** — colunas: ``label``, ``tipo`` (``SecaoBase``, ``Secao`` ou
    ``SecaoRetangular``), ``material`` (obrigatório para ``Secao`` e
    ``SecaoRetangular``), ``params`` (``chave=valor, chave=valor``).

    **nos** — colunas: ``id``, coordenadas e, opcionalmente, ``carga``,
    ``deslocamentos_nulos``, ``deslocamentos_prescritos``, ``apoio_elastico``
    (valores escritos como literais Python).

    **elementos** — colunas: ``id``, ``id_noI``, ``id_noJ``, ``sec`` e,
    opcionalmente, ``carga``, ``dTemp``, ``inclui_peso_proprio``.

    Parâmetros
    ----------
    caminho : str
        Caminho para o arquivo (``.xlsx``, ``.ods``, ``.xls``, …).
        Requer ``openpyxl`` para ``.xlsx`` e ``odfpy`` para ``.ods``.
    modelo : Modelo
        Instância do modelo a ser populada.
    """
    xls = pd.ExcelFile(caminho)

    for aba in ("materiais", "secoes", "nos", "elementos"):
        if aba not in xls.sheet_names:
            raise ValueError(
                f"Aba '{aba}' não encontrada no arquivo. "
                f"Abas disponíveis: {xls.sheet_names}"
            )

    materiais = ler_materiais(pd.read_excel(xls, sheet_name="materiais", dtype=str))
    secoes    = ler_secoes(pd.read_excel(xls, sheet_name="secoes", dtype=str), materiais)
    ler_nos(pd.read_excel(xls, sheet_name="nos"), modelo)
    ler_elementos(pd.read_excel(xls, sheet_name="elementos"), secoes, modelo)


# ---------------------------------------------------------------------------
# Escrita em planilha
# ---------------------------------------------------------------------------

def _coletar_materiais(modelo: Modelo) -> tuple[dict, list]:
    """Retorna ``(id_mat → label, [(label, Material)])`` para todas as seções com material."""
    from .secao import Secao

    vistos: dict[int, str] = {}   # id(mat) -> label
    ordem: list = []              # [(label, mat)]
    contador = 0
    for sec in modelo._secoes.values():
        if isinstance(sec, Secao) and id(sec.mat) not in vistos:
            contador += 1
            label = sec.mat.nome if sec.mat.nome else f"mat_{contador}"
            vistos[id(sec.mat)] = label
            ordem.append((label, sec.mat))
    return vistos, ordem


def _params_str(sec) -> str:
    """Converte os parâmetros de uma seção na string ``'chave=valor, ...'``."""
    from .secao import Secao, SecaoRetangular

    if isinstance(sec, SecaoRetangular):
        return f"b={sec.b!r}, h={sec.h!r}"
    if isinstance(sec, Secao):
        partes = [f"A={sec.A!r}"]
        if sec.I3: partes.append(f"I3={sec.I3!r}")
        if sec.I2: partes.append(f"I2={sec.I2!r}")
        if sec.J:  partes.append(f"J={sec.J!r}")
        return ", ".join(partes)
    # SecaoBase
    campos = [
        ("EA", sec.EA), ("EI3", sec.EI3), ("EI2", sec.EI2), ("GJ", sec.GJ),
        ("peso_unitario", sec.peso_unitario),
        ("cdtEA", sec.cdtEA), ("cdtEI2", sec.cdtEI2), ("cdtEI3", sec.cdtEI3),
    ]
    return ", ".join(f"{k}={v!r}" for k, v in campos if v)


def _df_materiais(modelo: Modelo, _id_mat: dict) -> pd.DataFrame:
    _, ordem = _coletar_materiais(modelo)
    linhas = []
    for label, mat in ordem:
        linhas.append({
            "label": label,
            "E":     mat.E,
            "G":     mat.G     if mat.G     else None,
            "pe":    mat.pe    if mat.pe    else None,
            "cdt":   mat.cdt   if mat.cdt   else None,
            "nome":  mat.nome  if mat.nome  else None,
        })
    return pd.DataFrame(linhas)


def _df_secoes(modelo: Modelo, id_mat: dict) -> pd.DataFrame:
    from .secao import Secao

    linhas = []
    for label, sec in modelo._secoes.items():
        tipo_nome = type(sec).__name__
        linha: dict = {
            "label":  label,
            "tipo":   tipo_nome,
            "params": _params_str(sec),
        }
        if isinstance(sec, Secao):
            linha["material"] = id_mat[id(sec.mat)]
        linhas.append(linha)

    df = pd.DataFrame(linhas)
    # garante ordem de colunas: label, tipo, material, params
    cols = [c for c in ("label", "tipo", "material", "params") if c in df.columns]
    return df[cols]


def _df_nos(modelo: Modelo) -> pd.DataFrame:
    eixos = modelo._tipos_no[modelo.tipo].eixos_globais
    linhas = []
    for id_no, no in modelo._nos.items():
        linha: dict = {"id": id_no}
        for i, eixo in enumerate(eixos):
            linha[eixo] = no.coor[i]
        if no.carga:
            linha["carga"] = str(no.carga)
        if no.deslocamentos_nulos:
            linha["deslocamentos_nulos"] = str(list(no.deslocamentos_nulos))
        if no.deslocamentos_prescritos:
            linha["deslocamentos_prescritos"] = str(no.deslocamentos_prescritos)
        if no.apoio_elastico:
            linha["apoio_elastico"] = str(no.apoio_elastico)
        linhas.append(linha)
    return pd.DataFrame(linhas)


def _df_elementos(modelo: Modelo) -> pd.DataFrame:
    no_para_id = {id(no): id_ext for id_ext, no in modelo._nos.items()}
    linhas = []
    for id_el, el in modelo._elementos.items():
        linha: dict = {
            "id":     id_el,
            "id_noI": no_para_id[id(el.noI)],
            "id_noJ": no_para_id[id(el.noJ)],
            "sec":    el.sec.nome,
        }
        if el.carga:
            linha["carga"] = str(el.carga)
        if el.dTemp:
            linha["dTemp"] = str(el.dTemp)
        if el.inclui_peso_proprio:
            linha["inclui_peso_proprio"] = True
        linhas.append(linha)
    return pd.DataFrame(linhas)


def _df_resultados_nos(modelo: Modelo) -> pd.DataFrame:
    from .graficos import conv_no

    assert modelo._res_nos is not None
    conv = modelo.conv
    linhas = []
    for id_no, res in modelo._res_nos.items():
        res_conv = conv_no(res, conv) if conv is not None else res
        linhas.append({"id": id_no, **res_conv})
    return pd.DataFrame(linhas)


def _df_resultados_elementos(modelo: Modelo) -> pd.DataFrame:
    from .graficos import conv_rel

    assert modelo._res_elementos is not None
    conv = modelo.conv
    linhas = []
    for id_el, res in modelo._res_elementos.items():
        res_convertido = conv_rel(res, conv) if conv is not None else res
        for grandeza, arr in res_convertido.items():
            linha: dict = {"elemento": id_el, "resultado": grandeza}
            for i, val in enumerate(arr):
                linha[f"pt_{i}"] = val
            linhas.append(linha)
    return pd.DataFrame(linhas)


def salvar_excel(caminho: str, modelo: Modelo) -> None:
    """Salva os dados de entrada e resultados do modelo em um arquivo Excel.

    Gera as seguintes abas:

    **materiais**, **secoes**, **nos**, **elementos** — dados de entrada no
    mesmo formato lido por :func:`carregar_excel`, permitindo recarregar o
    modelo a partir do arquivo gerado.

    **resultados_nos** — deslocamentos e reações em cada nó (apenas se o
    modelo foi analisado).

    **resultados_elementos** — esforços e deslocamentos ao longo de todos os
    elementos em uma única aba. A coluna ``elemento`` repete o ID do elemento
    para cada ponto; as demais colunas contêm as grandezas calculadas (``x``,
    ``N``, ``V2``, ``M3``, etc.).

    Parâmetros
    ----------
    caminho : str
        Caminho de destino (``.xlsx``). Requer ``openpyxl``.
    modelo : Modelo
        Instância do modelo (deve estar analisado para gerar as abas de
        resultados).
    """
    id_mat, _ = _coletar_materiais(modelo)

    with pd.ExcelWriter(caminho, engine="openpyxl") as writer:
        _df_materiais(modelo, id_mat).to_excel(writer, sheet_name="materiais", index=False)
        _df_secoes(modelo, id_mat).to_excel(writer, sheet_name="secoes",    index=False)
        _df_nos(modelo).to_excel(writer,            sheet_name="nos",       index=False)
        _df_elementos(modelo).to_excel(writer,      sheet_name="elementos", index=False)

        if modelo._res_nos is not None:
            _df_resultados_nos(modelo).to_excel(
                writer, sheet_name="resultados_nos", index=False
            )
        if modelo._res_elementos is not None:
            _df_resultados_elementos(modelo).to_excel(
                writer, sheet_name="resultados_elementos", index=False
            )
