"""
Módulo para definição dos elementos na análise estrutural
pelo Método da Rigidez Direta (MRD) no AECPy
"""

import numpy as np

from . import axial, flexao
from . import procedimentos as pmm
from .no import NoBase, NoPP, NoPE, NoTP, NoTE, NoGR
from .secao import Secao


def _acumular(dest, src):
    """Acumula os valores de ``src`` em ``dest``, somando os já existentes."""
    for key, value in src.items():
        if key in dest:
            dest[key] += value
        else:
            dest[key] = value


class ElementoBase:
    """Classe base de elemento para modelos estruturais via Método da Rigidez Direta.

    Cada subclasse define o tipo de modelo (PP, PE, TP, TE, GR), incluindo
    as partes de rigidez, cargas e temperatura disponíveis.
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
        self.__def_ini = None
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
        """Cargas distribuídas por direção, como ``(wi, wf)`` no início e no fim.

        Exemplo: ``{'w1': (10.0, 5.0), 'wz': -2.0}``
        """
        return dict(self.__carga)

    @carga.setter
    def carga(self, valores):
        """Define cargas distribuídas por dicionário; use None para limpar.

        Exemplo: ``{'w1': (10.0, 5.0), 'wz': -2.0}``
        """
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
        """Define se o peso próprio é incluído no carregamento do elemento."""
        if not isinstance(valor, bool):
            raise TypeError("inclui_peso_proprio deve ser do tipo bool")
        self.__inclui_peso_proprio = valor

    # ------------------------------------------------------------------
    # Variação de temperatura

    @property
    def dTemp(self):
        """Variações de temperatura por componente.

        - ``dT0``: variação uniforme de temperatura no centroide da seção.
        - ``dT2``: gradiente de temperatura na direção do eixo 2 (curvatura no plano 1-2).
        - ``dT3``: gradiente de temperatura na direção do eixo 3 (curvatura no plano 1-3).

        Exemplo: ``{'dT0': 20.0, 'dT2': -5.0}``
        """
        return dict(self.__dTemp)

    @dTemp.setter
    def dTemp(self, valores):
        """Define variações de temperatura por dicionário; use None para limpar.

        Exemplo: ``{'dT0': 20.0, 'dT2': -5.0}``
        """
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


    @property
    def def_ini(self):
        """Deformação inicial da barra, como um valor único para toda a extensão."""
        return self.__def_ini
    @def_ini.setter
    def def_ini(self, valor):
        """Define a deformação inicial da barra."""
        if not isinstance(valor, (int, float)):
            raise TypeError("def_ini deve ser um valor numérico")
        self.__def_ini = valor
    # ------------------------------------------------------------------
    # Cálculos

    def atualizar_geometria(self):
        """Recalcula o comprimento ``L`` e o vetor unitário ``e1`` de noI para noJ."""
        self.__L, self.__e1 = pmm.calc_L_u(self.noI.coor, self.noJ.coor)

    def R(self):
        """Matriz de transformação das coordenadas globais para as locais."""
        R = pmm.R3D(self.e1_3D)
        if self._R_ord is not None:
            R = pmm.reordenar_array(R, self._R_ord[0], self._R_ord[1])
        return R

    def Ke_local(self):
        """Matriz de rigidez nas coordenadas locais do elemento."""
        Kel = np.zeros((self.ngdl, self.ngdl))
        if "ax" in self._partes:
            pmm.espalhar(axial.calc_Ka(self.sec, self.L), Kel, self._partes["ax"])
        if "tr" in self._partes:
            pmm.espalhar(axial.calc_Kt(self.sec, self.L), Kel, self._partes["tr"])
        if "b3" in self._partes:
            pmm.espalhar(flexao.calc_Kb3(self.sec, self.L), Kel, self._partes["b3"])
        if "b2" in self._partes:
            pmm.espalhar(flexao.calc_Kb2(self.sec, self.L), Kel, self._partes["b2"])
        return Kel

    def carga_total(self) -> dict[str, np.ndarray]:
        """Cargas totais nos eixos locais do elemento.

        Converte cargas globais (wx, wy, wz) para locais (w1, w2, w3)
        e soma às cargas locais já definidas.
        """
        ct: dict[str, np.ndarray] = {c: np.array(self.__carga.get(c, (0.0, 0.0)))
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
        """Reações de engastamento perfeito devidas às cargas no elemento."""
        rep = np.zeros(self.ngdl)

        # Treliças (TP e TE): sem cargas distribuídas locais; peso próprio como
        # forças nodais equivalentes projetadas via matriz de rotação 3D
        if self.tipo in ("TP", "TE"):
            if self.inclui_peso_proprio:
                R = pmm.R3D(self.e1_3D)
                ff = R @ np.array([0.0, 0.0, self.sec.peso_unitario * self.L / 2])
                rep += np.concatenate((ff[:self.ndim], ff[:self.ndim]))
            dT = self.dTemp
            if "dT0" in dT:
                pmm.espalhar(axial.rep_T(self.L, self.sec, dT["dT0"]), rep, self._partes["ax"])
            if self.def_ini is not None:
                pmm.espalhar(axial.rep_def_ini(self.L, self.sec, self.def_ini), rep, self._partes["ax"])
            return rep  # demais partes não existem em treliças

        # Pórticos e grelhas: cargas locais via funções de REP por parte
        ct = self.carga_total()
        dT = self.dTemp
        if "ax" in self._partes:        # barra axial (PP, PE)
            pmm.espalhar(axial.rep_w1(self.L, ct["w1"]), rep, self._partes["ax"])
            if "dT0" in dT:
                pmm.espalhar(axial.rep_T(self.L, self.sec, dT["dT0"]), rep, self._partes["ax"])
            if self.def_ini is not None:
                pmm.espalhar(axial.rep_def_ini(self.L, self.sec, self.def_ini), rep, self._partes["ax"])
        if "b3" in self._partes:        # viga no plano 1-2 (PP, PE, GR)
            pmm.espalhar(flexao.rep_w2(self.L, ct["w2"]), rep, self._partes["b3"])
            if "dT2" in dT:
                pmm.espalhar(flexao.rep_T2(self.L, self.sec, dT["dT2"]), rep, self._partes["b3"])
        if "b2" in self._partes:        # viga no plano 1-3 (PE)
            pmm.espalhar(flexao.rep_w3(self.L, ct["w3"]), rep, self._partes["b2"])
            if "dT3" in dT:
                pmm.espalhar(flexao.rep_T3(self.L, self.sec, dT["dT3"]), rep, self._partes["b2"])
        return rep

    def rel(self, dl, x):
        """Deslocamentos e esforços internos na posição ``x`` do elemento."""
        resultado = {"x": x}
        xi = x / self.L
        pmm.check_xi(xi)
        ct = self.carga_total()
        dT = self.dTemp

        if "ax" in self._partes:        # barra axial (PP, PE, TP, TE)
            il = self._partes["ax"]
            _acumular(resultado, axial.rel_d_a(xi, self.sec, self.L, dl[il]))
            if "w1" in ct:              # treliças não possuem carga local w1
                _acumular(resultado, axial.rel_w1(xi, self.sec, self.L, ct["w1"]))
            if "dT0" in dT:
                _acumular(resultado, axial.rel_T(xi, self.sec, self.L, dT["dT0"]))
            if self.def_ini is not None:
                _acumular(resultado, axial.rel_def_ini(xi, self.sec, self.L, self.def_ini))
        if "tr" in self._partes:        # torção (PE)
            il = self._partes["tr"]
            _acumular(resultado, axial.rel_d_t(xi, self.sec, self.L, dl[il]))
        if "b3" in self._partes:        # viga no plano 1-2 (PP, PE, GR)
            il = self._partes["b3"]
            _acumular(resultado, flexao.rel_d_b3(xi, self.sec, self.L, dl[il]))
            _acumular(resultado, flexao.rel_w2(xi, self.sec, self.L, ct["w2"]))
            if "dT2" in dT:
                _acumular(resultado, flexao.rel_T2(xi, self.sec, self.L, dT["dT2"]))
        if "b2" in self._partes:        # viga no plano 1-3 (PE)
            il = self._partes["b2"]
            _acumular(resultado, flexao.rel_d_b2(xi, self.sec, self.L, dl[il]))
            _acumular(resultado, flexao.rel_w3(xi, self.sec, self.L, ct["w3"]))
            if "dT3" in dT:
                _acumular(resultado, flexao.rel_T3(xi, self.sec, self.L, dT["dT3"]))
        return resultado

    def dl(self, d_global):
        """Deslocamentos nodais nos eixos locais a partir do vetor de deslocamentos global."""
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
        """Define a alteração dos eixos locais: ângulo, NoBase ou coordenadas de um ponto."""
        if isinstance(alt, (float, int, NoBase)):
            self.__alterar_el = alt
        elif len(alt) == self.noI.ndim:
            self.__alterar_el = np.array(alt)
        else:
            raise ValueError(
                "alt deve ser um ângulo (em graus), um NoBase ou coordenadas de um ponto"
            )

    def R(self):
        """Matriz de transformação considerando a alteração dos eixos locais."""
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
