"""
Módulo para definição dos elementos na análise estrutural
pelo Método da Rigidez Direta (MRD) no AECPy
"""

import numpy as np

from . import axial, flexao
from . import procedimentos as pmm
from .no import No
from .secao import Secao


class Elemento:
    """
    Elemeto para análise estrutural pelo Método da Rigidez Direta no AECPy

    Atributes
    ---------
    noI, noJ: No
        Nós inicial e final do elemento
    sec: Secao
        Seção transversal do elemento
    mat: Material
        Material que forma o elemento
    L: float
        Comprimento do elemento
    e1: numpy.ndarray (shape = (3,))
        Vetor unitário na direção do noI para o noJ
    """

    __parte_rigidez = {
        "ax": {
            "calc_K": axial.calc_Ka,
            "rel_d": axial.rel_d_a,
            "carga": "w1",
            "rep": axial.rep_w1,
            "rel_w": axial.rel_w1,
            "dTemp": "dT0",
            "rep_dT": axial.rep_T,
            "rel_dT": axial.rel_T,
        },
        "tr": {
            "calc_K": axial.calc_Kt,
            "rel_d": axial.rel_d_t,
            "carga": "",
            "rep": None,
            "rel_w": None,
            "dTemp": "",
            "rep_dT": None,
            "rel_dT": None,
        },
        "b2": {
            "calc_K": flexao.calc_Kb2,
            "rel_d": flexao.rel_d_b2,
            "carga": "w3",
            "rep": flexao.rep_w3,
            "rel_w": flexao.rel_w3,
            "dTemp": "dT3",
            "rep_dT": flexao.rep_T3,
            "rel_dT": flexao.rel_T3,
        },
        "b3": {
            "calc_K": flexao.calc_Kb3,
            "rel_d": flexao.rel_d_b3,
            "carga": "w2",
            "rep": flexao.rep_w2,
            "rel_w": flexao.rel_w2,
            "dTemp": "dT2",
            "rep_dT": flexao.rep_T2,
            "rel_dT": flexao.rel_T2,
        },
    }

    __carga = dict()
    __incluir_peso_proprio = False
    __dTemp = dict()

    # __quantidades = {   'força': ['N','V2','V3'] + [f'f{i}' for i in [1,2,3,'x','y','z']],
    #                     'momento':['Mt','M2','M3'] + [f'm{i}' for i in [1,2,3,'x','y','z']],
    #                     'comprimento':['x','y','z'],
    #                     'translação':[f'u{i}' for i in [1,2,3,'x','y','z']],
    #                     'rotação': [f'r{i}' for i in [1,2,3,'x','y','z']]
    #                     }
    # __unidade_padrao_saida = {  'força': 'kN',
    #                             'momento': 'kNm',
    #                             'comprimento': 'cm',
    #                             'translação': 'mm',
    #                             'rotação': 'rad'}

    @classmethod
    def iniciar(cls, tipo):
        if tipo == "TP":  # elemento de treliça plana
            cls.__R_ord = ([0, 1], (0, 2))  # (u1,u2) = R (ux,uz)
            cls.__partes = {"ax": [0, 2]}
            cls.__cargas_locais = tuple()
            cls.__cargas_globais = tuple()
            cls.__var_temperatura = ("dT0",)

        elif tipo == "TE":  # elemento de treliça espacial
            cls.__R_ord = None
            cls.__partes = {"ax": [0, 3]}
            cls.__cargas_locais = tuple()
            cls.__cargas_globais = tuple()
            cls.__var_temperatura = ("dT0",)

        elif tipo == "PP":  # elemento de pórtico plano
            cls.__R_ord = ([0, 1, 2], [0, 2, 1])  # (u1,u2,r3) = R (ux,uz,ry)
            cls.__partes = {"ax": [0, 3], "b3": [1, 2, 4, 5]}
            cls.__cargas_locais = ("w1", "w2")
            cls.__cargas_globais = ("wx", "wz")
            cls.__var_temperatura = ("dT0", "dT2")

        elif tipo == "GR":  # elemento de grelha
            cls.__R_ord = ([1, 0, 2], [2, 0, 1])  # (u2,r1,r3) = R (uz,rx,ry)
            cls.__partes = {"tr": [1, 4], "b3": [0, 2, 3, 5]}
            cls.__cargas_locais = ("w2",)
            cls.__cargas_globais = ("wz",)
            cls.__var_temperatura = ("dT2",)

        elif tipo == "PE":  # elemento de pórtico espacial
            cls.__R_ord = None
            cls.__partes = {
                "ax": [0, 6],
                "tr": [3, 9],
                "b3": [1, 5, 7, 11],
                "b2": [2, 4, 8, 10],
            }
            cls.__cargas_locais = ("w1", "w2", "w3")
            cls.__cargas_globais = ("wx", "wy", "wz")
            cls.__var_temperatura = ("dT0", "dT2", "dT3")

        cls.__tipo = tipo
        cls.__cargas_permitidas = cls.__cargas_globais + cls.__cargas_locais

    @classmethod
    def incluir_peso_proprio(cls, incluir):
        """Define a inclusão do peso próprio no carregamento do elemento"""
        if isinstance(incluir, bool):
            cls.__incluir_peso_proprio = incluir
        else:
            raise ValueError("incluir deve ser bool")

    def __init__(self, noI, noJ, sec) -> None:
        if (not isinstance(noI, No)) or (not isinstance(noJ, No)):
            raise TypeError("noI ou noJ não são do tipo `No`")
        if not isinstance(sec, Secao):
            raise TypeError("sec não é do tipo `Secao`")
        self.noI = noI  # nó inicial
        self.noJ = noJ  # nó final
        self.sec = sec  # seçaõ transvreal
        self.mat = sec.mat  # matreal do elemento

        # cálculo de comprimento e direção do elemento
        self.atualizar_geometria()

    @property
    def tipo(self):
        return self.__tipo

    @property
    def ndim(self):
        """Número de dimensões do espaço no modelo (2 ou 3)"""

        return self.noI.ndim

    @property
    def ngdl(self):
        """Número de gdl do elemento"""

        return 2 * self.noI.ngdl

    @property
    def igdl(self):
        """Índices globais dos gdl do elemento"""

        return list(self.noI.igdl) + list(self.noJ.igdl)

    @property
    def L(self):
        """Comprimento do elemento"""
        return self.__L

    @property
    def e1(self):
        """Vetor unitário na direção do noI para noJ"""
        return self.__e1

    @property
    def e1_3D(self):
        """Vetor unitário na direção do noI para noJ, no espaço 3D"""
        if self.ndim == 2:
            if self.tipo in ["TP", "PP"]:
                return np.array([self.e1[0], 0, self.e1[1]])
            # self.tipo == 'GR'
            return np.array([self.e1[0], self.e1[1], 0])
        return self.e1

    @property
    def carga(self):
        """
        Define os carregamentos no elemento
        Retorna um dicionário onde
            > key: direção de cada carregamento
            impo> value: (wi,wf) com a intensidade da carga distribuída no início e no fim do elemento"""
        return self.__carga

    @property
    def carregado(self):
        return len(self.__carga) > 0

    @property
    def dTemp(self):
        return self.__dTemp

    @property
    def variacao_termica(self):
        return len(self.__dTemp) > 0

    @property
    def alterar_eixos_locais(self):
        """
        Dado para alteração dos eixos locais
        se float ou int: representa o ângulo de rotação dos eixos 2,3 em graus e no sentido anti-horário
        se No ou array: representa um ponto adicional usado para criar um plano que contém os eixos 1-2
        """
        try:
            return self.__alterar_el
        except AttributeError:
            return None  # if hasattr(self, '_Elemento_alterar_el') else None

    @alterar_eixos_locais.setter
    def alterar_eixos_locais(self, alt):
        if self.tipo != "PE":
            raise ValueError(
                "A alteração dos eixos locais só é possível em Pórticos Espaciais"
            )
        if isinstance(alt, (float, int, No)):
            self.__alterar_el = alt
        elif len(alt) == self.noI.ndim:
            self.__alterar_el = np.array(alt)
        else:
            raise ValueError(
                "alt deve ser um ângulo (em rad), um No ou coordenadas de um ponto"
            )

    def atualizar_geometria(self):
        """Calcula o comprimento (L) e o vetor unitário de noI par noJ (e1)"""

        self.__L, self.__e1 = pmm.calc_L_u(self.noI.coor, self.noJ.coor)

    def R(self):
        """
        Calcula a matriz de transformação das coordenadas globais para as locais do elemento
        """
        alt = self.alterar_eixos_locais
        if alt:
            # Eixos locais em direção diferente da padrão
            if isinstance(alt, float):
                # rotação dos eixos 2-3 definida por um ângulo
                Rpadrao = pmm.R3D(self.e1_3D)
                R = pmm.R3D_mod_ang(Rpadrao, alt)
            else:
                # direção do eixo 2 definda no plano e1 , u
                if isinstance(alt, No):
                    # vetor u definido por um 3º nó
                    u = alt.coor - self.noI.coor
                else:
                    # vetor u definido por um ponto adicional
                    u = alt - self.noI.coor
                u /= np.linalg.norm(u)
                R = pmm.R3D_u(self.e1_3D, u)
        else:
            # eixos locais na direção padrão
            R = pmm.R3D(self.e1_3D)

        if self.__R_ord is not None:
            R = pmm.reordenar_array(R, self.__R_ord[0], self.__R_ord[1])

        return R

    def Ke_local(self):
        """Calcula as parcelas de rigidez do elemento e monta a matriz de rigidz local"""
        ngdl = self.ngdl
        Kel = np.zeros((ngdl, ngdl))
        for prt, il in self.__partes.items():
            calc_K = self.__parte_rigidez[prt]["calc_K"]
            pmm.espalhar(calc_K(self.sec, self.L), Kel, il)
        return Kel

    def definir_carga(self, **kwargs):
        """Introdução das cargas no elemento"""
        # forças extrenas distribuídas em todo o elemento
        self.__carga = dict()
        for carga, w in kwargs.items():
            if carga in self.__cargas_permitidas:
                if isinstance(w, (int, float)):  # carregamento constante
                    self.__carga[carga] = (w, w)
                elif len(w) == 2 and all(
                    [isinstance(v, (int, float)) for v in w]
                ):
                    self.__carga[carga] = tuple(w)  # carregamento linear
                else:
                    raise ValueError(f"A carga {carga} tem formato incorreto")
            else:
                raise ValueError(
                    f"A carga {carga} não está entra as permitidas: {self.__cargas_permitidas}"
                )
        return

    def definir_carga_termica(self, **kwargs):
        """Introdução da variação da temperatura no elemento"""
        # variação de temperatura em todo o elemento
        self.__dTemp = dict()
        for dTemp, T in kwargs.items():
            if dTemp in self.__var_temperatura:
                if isinstance(T, (int, float)):
                    # variação de termperatura uniforme na seção transversal
                    self.__dTemp[dTemp] = T

                else:
                    raise ValueError(
                        f"A variação de tempertaura {dTemp} tem formato incorreto: {T}"
                    )
            else:
                raise ValueError(
                    f"A variação de temperatura {dTemp} não está entre as perimtidas: {self.__var_temperatura}"
                )

        return

    def carga_total(self):
        """
        Cálculo das cargas totais no elemento, somando as cargas
        Transforma cargas dadas segundo coordenadas globais (wx,wy,wz)
        em cagas segundo coordenadas locais (w1,w2,w3) e soma às cargas
        dadas originalmente segundo os eixos locais.
        """
        carga_total = dict()

        # cargas locais
        for c in self.__cargas_locais:
            if c in self.carga:  # carga local não nula
                carga_total[c] = np.array(self.carga[c])
            else:  # carga local nula
                carga_total[c] = np.zeros(2)

        # Conversão de cargas globais em locais

        wg_i = np.zeros(3)  # intensidade no início do elemento [wx,wy,wz]
        wg_j = np.zeros(3)  # intensidade no fim do elemento
        for i, c in enumerate(["wx", "wy", "wz"]):
            if c in self.carga:
                wg_i[i] = self.carga[c][0]
                wg_j[i] = self.carga[c][1]

        if self.__incluir_peso_proprio:
            wg_i[2] -= self.sec.peso_unitario
            wg_j[2] -= self.sec.peso_unitario

        # Transformação das coordenadas globais para as locais
        R = pmm.R3D(self.e1_3D)
        wl_i = R @ wg_i
        wl_j = R @ wg_j

        for i, c in enumerate(["w1", "w2", "w3"]):
            if c in self.__cargas_locais:
                carga_total[c][0] += wl_i[i]
                carga_total[c][1] += wl_j[i]

        return carga_total


    def rep(self):
        """Calcula as reações de engastamento perfeito das cargas nos elementos"""

        rep = np.zeros(self.ngdl)
        carga_total = self.carga_total()

        # <<<<< atenação --- programar para treliças ---- >>>>>
        for prt, il in self.__partes.items():
            carga = self.__parte_rigidez[prt]["carga"]
            calc_rep = self.__parte_rigidez[prt]["rep"]
            if calc_rep is not None:
                w = carga_total[carga]
                rep[il] += calc_rep(self.L, w)

            dTemp = self.__parte_rigidez[prt]["dTemp"]
            calc_rep_dT = self.__parte_rigidez[prt]["rep_dT"]
            if calc_rep_dT is not None and dTemp in self.__dTemp:
                dT = self.__dTemp[dTemp]
                rep[il] += calc_rep_dT(self.L, self.sec, dT)

        return rep

    def rel(self, dl, x):
        """Calcula os deslocamentos e esforços internos em uma posição x do elemento"""

        rel = {"x": x}

        xi = x / self.L
        pmm.check_xi(xi)

        carga_total = self.carga_total()

        # esforços e deslocamentos no elemento
        for prt, il in self.__partes.items():
            # devido aos deslocamentos nodais
            dl_prt = dl[il]
            calc_rel_d = self.__parte_rigidez[prt]["rel_d"]
            rel_d = calc_rel_d(xi, self.sec, self.L, dl_prt)

            # esforços e deslocamentos no elemento devido aos carregamentos externos
            calc_rel_w = self.__parte_rigidez[prt]["rel_w"]
            if calc_rel_w is not None:
                carga = self.__parte_rigidez[prt]["carga"]
                w_prt = carga_total.get(carga,[0.0,0.0])
                rel_w = calc_rel_w(xi, self.sec, self.L, w_prt)
            else:
                rel_w = None

            # esforços e deslocamentos no elemento devido à variação de temperatura
            calc_rel_dT = self.__parte_rigidez[prt]["rel_dT"]
            dTemp = self.__parte_rigidez[prt]["dTemp"]
            if calc_rel_dT is not None and dTemp in self.__dTemp:
                dT_prt = self.__dTemp[dTemp]
                rel_dT = calc_rel_dT(xi, self.sec, self.L, dT_prt)
            else:
                rel_dT = None
            # somando os rel da parte no rel do elemento
            for dd in [rel_d, rel_w, rel_dT]:
                if dd is None:
                    continue
                for key, value in dd.items():
                    if key in rel:
                        rel[key] += value
                    else:
                        rel[key] = value
        return rel

    def dl(self, d_global):
        """
        Calcula os deslocamentos dos nós do elemento em coordenadas
        locais a partir do vetor global de deslocamentos em coordenadas
        globais:

        Parâmetros:
        ----------
        d_global: np.ndarray
            vetor global de deslocamentos em coordenadas globais
        """
        dg = d_global[self.igdl]
        R = self.R()
        return pmm.transf_coord(dg, R)
