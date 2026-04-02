"""
Módulo para definição dos elementos na análise estrutural
pelo Método da Rigidez Direta (MRD) no AECPy
"""

import numpy as np

from . import axial, flexao
from . import procedimentos as pmm
from .no import NoBase, NoPP, NoPE, NoTP, NoTE, NoGR
from .secao import Secao


class ElementoBase:
    """Classe base de elemento para modelos estruturais via Método da Rigidez Direta.

    Cada subclasse define o tipo de modelo (PP, PE, TP, TE, GR), incluindo
    as partes de rigidez, cargas e temperatura disponíveis.

    Atributos
    ---------
    noI, noJ : NoBase
        Nós inicial e final do elemento.
    sec : Secao
        Seção transversal do elemento.
    L : float
        Comprimento do elemento.
    e1 : numpy.ndarray (shape = (3,))
        Vetor unitário na direção do noI para o noJ.
    """

    # atributos do tipo estrutural (definidos em subclasses)
    tipo = ""
    ndim = 0                # dimensão do espaço (2 ou 3)
    ngdl = 0                # número de gdl do elemento (2 × ngdl do nó)
    _tipo_no = None         # classe de nó compatível com o elemento
    _R_ord = None           # reordenamento da matriz R; None para modelos 3D
    _partes = {}            # partes de rigidez e seus índices de gdl locais
    _cargas_locais = ()     # identificadores das cargas nos eixos locais
    _cargas_globais = ()    # identificadores das cargas nos eixos globais
    _var_temperatura = ()   # identificadores das variáveis de temperatura

    # mapa de funções de rigidez por parte (compartilhado entre todos os tipos)
    _parte_rigidez = {
        "ax": {
            "calc_K": axial.calc_Ka,
            "rel_d":  axial.rel_d_a,
            "carga":  "w1",
            "rep":    axial.rep_w1,
            "rel_w":  axial.rel_w1,
            "dTemp":  "dT0",
            "rep_dT": axial.rep_T,
            "rel_dT": axial.rel_T,
        },
        "tr": {
            "calc_K": axial.calc_Kt,
            "rel_d":  axial.rel_d_t,
            "carga":  "",
            "rep":    None,
            "rel_w":  None,
            "dTemp":  "",
            "rep_dT": None,
            "rel_dT": None,
        },
        "b2": {
            "calc_K": flexao.calc_Kb2,
            "rel_d":  flexao.rel_d_b2,
            "carga":  "w3",
            "rep":    flexao.rep_w3,
            "rel_w":  flexao.rel_w3,
            "dTemp":  "dT3",
            "rep_dT": flexao.rep_T3,
            "rel_dT": flexao.rel_T3,
        },
        "b3": {
            "calc_K": flexao.calc_Kb3,
            "rel_d":  flexao.rel_d_b3,
            "carga":  "w2",
            "rep":    flexao.rep_w2,
            "rel_w":  flexao.rel_w2,
            "dTemp":  "dT2",
            "rep_dT": flexao.rep_T2,
            "rel_dT": flexao.rel_T2,
        },
    }

    def __init__(self, noI, noJ, sec) -> None:
        if self.tipo == "":
            raise TypeError(
                "ElementoBase é base. Use uma subclasse como "
                "ElementoPE, ElementoPP, ElementoTE, ElementoTP ou ElementoGR."
            )
        assert self._tipo_no is not None  # garantido pelas subclasses
        if not isinstance(noI, self._tipo_no) or not isinstance(noJ, self._tipo_no):
            raise TypeError(
                f"noI e noJ devem ser do tipo {self._tipo_no.__name__} "
                f"para o elemento do tipo '{self.tipo}'"
            )
        if not isinstance(sec, Secao):
            raise TypeError("sec não é do tipo Secao")

        self.noI = noI
        self.noJ = noJ
        self.sec = sec

        self.__carga = {}
        self.__dTemp = {}
        self.__inclui_peso_proprio = False

        self.atualizar_geometria()

    # ------------------------------------------------------------------
    # Propriedades geométricas e estruturais

    @property
    def igdl(self):
        """Índices globais dos gdl do elemento."""
        return list(self.noI.igdl) + list(self.noJ.igdl)

    @property
    def L(self):
        """Comprimento do elemento."""
        return self.__L

    @property
    def e1(self):
        """Vetor unitário na direção do noI para o noJ."""
        return self.__e1

    @property
    def e1_3D(self):
        """Vetor unitário na direção do noI para o noJ, no espaço 3D."""
        if self.ndim == 2:
            if self.tipo in ("TP", "PP"):
                return np.array([self.e1[0], 0.0, self.e1[1]])
            # GR: plano xOy
            return np.array([self.e1[0], self.e1[1], 0.0])
        return self.e1

    # ------------------------------------------------------------------
    # Carregamento

    @property
    def carga(self):
        """Dicionário de carregamentos distribuídos no elemento.

        Cada entrada associa a direção à intensidade ``(wi, wf)`` no início
        e no fim do elemento. Atribuir ``None`` limpa os carregamentos.

        Cargas permitidas dependem do tipo: eixos globais (ex.: ``wx``, ``wz``)
        e locais (ex.: ``w1``, ``w2``).
        """
        return dict(self.__carga)

    @carga.setter
    def carga(self, valores):
        if valores is None:
            self.__carga = {}
            return
        if not isinstance(valores, dict):
            raise TypeError("carga deve ser um dicionário ou None")
        cargas_permitidas = self._cargas_globais + self._cargas_locais
        nova_carga = {}
        for nome, w in valores.items():
            if nome not in cargas_permitidas:
                raise ValueError(
                    f"A carga '{nome}' não está entre as permitidas: {cargas_permitidas}"
                )
            if isinstance(w, (int, float)):
                nova_carga[nome] = (w, w)
            elif len(w) == 2 and all(isinstance(v, (int, float)) for v in w):
                nova_carga[nome] = tuple(w)
            else:
                raise ValueError(f"A carga '{nome}' tem formato incorreto")
        self.__carga = nova_carga

    @property
    def inclui_peso_proprio(self):
        """``True`` se o peso próprio é incluído no carregamento do elemento."""
        return self.__inclui_peso_proprio

    @inclui_peso_proprio.setter
    def inclui_peso_proprio(self, valor):
        if not isinstance(valor, bool):
            raise TypeError("inclui_peso_proprio deve ser do tipo bool")
        self.__inclui_peso_proprio = valor

    # ------------------------------------------------------------------
    # Variação de temperatura

    @property
    def dTemp(self):
        """Dicionário de variações de temperatura no elemento.

        Atribuir ``None`` limpa as variações. Variáveis permitidas dependem
        do tipo de elemento (ex.: ``dT0``, ``dT2``, ``dT3``).
        """
        return dict(self.__dTemp)

    @dTemp.setter
    def dTemp(self, valores):
        if valores is None:
            self.__dTemp = {}
            return
        if not isinstance(valores, dict):
            raise TypeError("dTemp deve ser um dicionário ou None")
        for var, T in valores.items():
            if var not in self._var_temperatura:
                raise ValueError(
                    f"A variação de temperatura '{var}' não está entre as "
                    f"permitidas: {self._var_temperatura}"
                )
            if not isinstance(T, (int, float)):
                raise TypeError(f"O valor de temperatura '{var}' deve ser numérico")
        self.__dTemp = dict(valores)

    # ------------------------------------------------------------------
    # Cálculos

    def atualizar_geometria(self):
        """Calcula o comprimento (L) e o vetor unitário de noI para noJ (e1)."""
        self.__L, self.__e1 = pmm.calc_L_u(self.noI.coor, self.noJ.coor)

    def R(self):
        """Calcula a matriz de transformação das coordenadas globais para as locais."""
        R = pmm.R3D(self.e1_3D)
        if self._R_ord is not None:
            R = pmm.reordenar_array(R, self._R_ord[0], self._R_ord[1])
        return R

    def Ke_local(self):
        """Calcula a matriz de rigidez local do elemento."""
        Kel = np.zeros((self.ngdl, self.ngdl))
        for prt, il in self._partes.items():
            calc_K = self._parte_rigidez[prt]["calc_K"]
            pmm.espalhar(calc_K(self.sec, self.L), Kel, il)
        return Kel

    def carga_total(self):
        """Cargas totais em coordenadas locais.

        Converte cargas em coordenadas globais (wx, wy, wz) para locais
        (w1, w2, w3) e soma às cargas definidas originalmente nos eixos locais.
        """
        ct = {c: np.array(self.__carga.get(c, (0.0, 0.0)))
              for c in self._cargas_locais}

        wg_i = np.zeros(3)
        wg_j = np.zeros(3)
        for i, c in enumerate(("wx", "wy", "wz")):
            if c in self.__carga:
                wg_i[i] = self.__carga[c][0]
                wg_j[i] = self.__carga[c][1]

        if self.__inclui_peso_proprio:
            wg_i[2] -= self.sec.peso_unitario
            wg_j[2] -= self.sec.peso_unitario

        R = pmm.R3D(self.e1_3D)
        wl_i = R @ wg_i
        wl_j = R @ wg_j

        for i, c in enumerate(("w1", "w2", "w3")):
            if c in self._cargas_locais:
                ct[c][0] += wl_i[i]
                ct[c][1] += wl_j[i]

        return ct

    def rep(self):
        """Calcula as reações de engastamento perfeito das cargas no elemento."""
        rep = np.zeros(self.ngdl)

        # caso especial: treliças com peso próprio
        if self.tipo.startswith("T") and self.__inclui_peso_proprio:
            R = pmm.R3D(self.e1_3D)
            ff = R @ np.array([0.0, 0.0, self.sec.peso_unitario * self.L / 2])
            if self.ndim == 2:
                return np.concatenate((ff[:2], ff[:2]))
            return np.concatenate((ff, ff))

        ct = self.carga_total()

        for prt, il in self._partes.items():
            calc_rep = self._parte_rigidez[prt]["rep"]
            if calc_rep is not None:
                carga = self._parte_rigidez[prt]["carga"]
                rep[il] += calc_rep(self.L, ct[carga])

            calc_rep_dT = self._parte_rigidez[prt]["rep_dT"]
            dTemp_key = self._parte_rigidez[prt]["dTemp"]
            if calc_rep_dT is not None and dTemp_key in self.__dTemp:
                rep[il] += calc_rep_dT(self.L, self.sec, self.__dTemp[dTemp_key])

        return rep

    def rel(self, dl, x):
        """Calcula deslocamentos e esforços internos na posição ``x`` do elemento."""
        rel = {"x": x}
        xi = x / self.L
        pmm.check_xi(xi)

        ct = self.carga_total()

        for prt, il in self._partes.items():
            rel_d = self._parte_rigidez[prt]["rel_d"](xi, self.sec, self.L, dl[il])

            calc_rel_w = self._parte_rigidez[prt]["rel_w"]
            if calc_rel_w is not None:
                carga = self._parte_rigidez[prt]["carga"]
                rel_w = calc_rel_w(xi, self.sec, self.L, ct.get(carga, [0.0, 0.0]))
            else:
                rel_w = None

            calc_rel_dT = self._parte_rigidez[prt]["rel_dT"]
            dTemp_key = self._parte_rigidez[prt]["dTemp"]
            if calc_rel_dT is not None and dTemp_key in self.__dTemp:
                rel_dT = calc_rel_dT(xi, self.sec, self.L, self.__dTemp[dTemp_key])
            else:
                rel_dT = None

            for dd in (rel_d, rel_w, rel_dT):
                if dd is None:
                    continue
                for key, value in dd.items():
                    if key in rel:
                        rel[key] += value
                    else:
                        rel[key] = value

        return rel

    def dl(self, d_global):
        """Deslocamentos nodais em coordenadas locais a partir do vetor global."""
        return pmm.transf_coord(d_global[self.igdl], self.R())


# ----------------------------------------------------------------------
# Subclasses por tipo estrutural
# ----------------------------------------------------------------------

class ElementoPP(ElementoBase):
    """Elemento para pórtico plano."""

    tipo = "PP"
    ndim = 2
    ngdl = 6
    _tipo_no = NoPP
    _R_ord = ([0, 1, 2], [0, 2, 1])        # (u1, u2, r3) = R (ux, uz, ry)
    _partes = {"ax": [0, 3], "b3": [1, 2, 4, 5]}
    _cargas_locais = ("w1", "w2")
    _cargas_globais = ("wx", "wz")
    _var_temperatura = ("dT0", "dT2")


class ElementoPE(ElementoBase):
    """Elemento para pórtico espacial."""

    tipo = "PE"
    ndim = 3
    ngdl = 12
    _tipo_no = NoPE
    _R_ord = None
    _partes = {"ax": [0, 6], "tr": [3, 9], "b3": [1, 5, 7, 11], "b2": [2, 4, 8, 10]}
    _cargas_locais = ("w1", "w2", "w3")
    _cargas_globais = ("wx", "wy", "wz")
    _var_temperatura = ("dT0", "dT2", "dT3")

    @property
    def alterar_eixos_locais(self):
        """Dado para alteração dos eixos locais.

        - ``float`` ou ``int``: ângulo de rotação dos eixos 2-3 em graus,
          sentido anti-horário.
        - ``NoBase`` ou array: ponto adicional para definir o plano
          que contém os eixos 1-2.
        """
        try:
            return self.__alterar_el
        except AttributeError:
            return None

    @alterar_eixos_locais.setter
    def alterar_eixos_locais(self, alt):
        if isinstance(alt, (float, int, NoBase)):
            self.__alterar_el = alt
        elif len(alt) == self.noI.ndim:
            self.__alterar_el = np.array(alt)
        else:
            raise ValueError(
                "alt deve ser um ângulo (em graus), um NoBase ou coordenadas de um ponto"
            )

    def R(self):
        """Calcula a matriz de transformação das coordenadas globais para as locais."""
        alt = self.alterar_eixos_locais
        if alt is None:
            return pmm.R3D(self.e1_3D)
        if isinstance(alt, (int, float)):
            return pmm.R3D_mod_ang(pmm.R3D(self.e1_3D), alt)
        u = (alt.coor if isinstance(alt, NoBase) else alt) - self.noI.coor
        return pmm.R3D_u(self.e1_3D, u / np.linalg.norm(u))


class ElementoTP(ElementoBase):
    """Elemento para treliça plana."""

    tipo = "TP"
    ndim = 2
    ngdl = 4
    _tipo_no = NoTP
    _R_ord = ([0, 1], (0, 2))               # (u1, u2) = R (ux, uz)
    _partes = {"ax": [0, 2]}
    _cargas_locais = ()
    _cargas_globais = ()
    _var_temperatura = ("dT0",)


class ElementoTE(ElementoBase):
    """Elemento para treliça espacial."""

    tipo = "TE"
    ndim = 3
    ngdl = 6
    _tipo_no = NoTE
    _R_ord = None
    _partes = {"ax": [0, 3]}
    _cargas_locais = ()
    _cargas_globais = ()
    _var_temperatura = ("dT0",)


class ElementoGR(ElementoBase):
    """Elemento para grelha."""

    tipo = "GR"
    ndim = 2
    ngdl = 6
    _tipo_no = NoGR
    _R_ord = ([1, 0, 2], [2, 0, 1])        # (u2, r1, r3) = R (uz, rx, ry)
    _partes = {"tr": [1, 4], "b3": [0, 2, 3, 5]}
    _cargas_locais = ("w2",)
    _cargas_globais = ("wz",)
    _var_temperatura = ("dT2",)
