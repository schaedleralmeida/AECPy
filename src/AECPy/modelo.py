"""
Módulo para definição do modelo da estrutura para análise
pelo Método da Rigidez Direta (MRD) no AECPy
"""

import numpy as np
import pandas as pd

from . import procedimentos as pmm
from .elemento import Elemento
from .no import No
from .unidades import Conversor
from .graficos import unidade_padrao_saida


def numerar_gdls(nos: list[No]) -> int:
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
        iS = sorted(list(no.il_gdlr) + list(no.il_gdlp))
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
    nos: list[No], els: list[Elemento], ngdl: int
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
        # if el.carregado or el.variacao_termica or el.def_ini:
        # Reações de engastamento perfeito (rep) em coordenadas locais
        rep_el_local = el.rep()
        # rep em coordenadas globais
        rep_el = pmm.transf_coord(rep_el_local, R, inv=True)
        # somando ao vetor de forças nodais equivalentes
        F[igdl] += -rep_el

    # Contribuição dos nós à matriz de rigidez global e vetor de forças nodais
    for no in nos:
        # nó com força externa
        if no.carregado:
            F[no.igdl] += no.p
        # nó com deslocamento prescrito
        if no.ngdlp > 0:
            for ig, dp in zip(no.igdlp, no.dp):
                F -= K[:, ig] * dp
        # nó com apoio elástico
        if no.ngdle > 0:
            for ig, k in zip(no.igdle, no.k):
                K[ig, ig] += k

    return K, F


def resolver_SEL(K: np.ndarray, F: np.ndarray, nos: list[No], ngdlF: int) -> tuple[
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
        if no.ngdlp > 0:
            U[no.igdlp] = no.dp

    # Reações
    R[ngdlF:] = K[ngdlF:, :ngdlF] @ U[:ngdlF] - F[ngdlF:]

    # Apoios elásticos
    for no in nos:
        if no.ngdle > 0:
            for ig, k_ae in zip(no.igdle, no.k):
                R[ig] += -k_ae * U[ig]

    return U, R


def resultados_nos(nos: list[No], U: np.ndarray, R: np.ndarray) -> pd.DataFrame:
    """
    Função para criar uma dataframe com os resultados nos nós
    """
    deslocamentos = [d for d in No.info("gdls globais")]
    reacoes = ["R" + f for f in No.info("forças globais")]
    ngdl_no = nos[0].ngdl
        
    data = np.zeros((len(nos), 2 * ngdl_no))

    for i, no in enumerate(nos):
        igdl = no.igdl
        data[i, :ngdl_no] = U[igdl]
        data[i, ngdl_no:] = R[igdl]

    df = pd.DataFrame(data, columns=(deslocamentos + reacoes))

    return df


def resultados_elementos(els: list[Elemento], U: np.ndarray, npts:int=5) -> list[dict]:
    """
    Função para criar uma dataframe com os resultados nos elementos
    """
    resultados = list()
    for el in els:
        x = np.linspace(0,el.L,npts)
        dl = el.dl(U)
        rel = el.rel(dl, x)
        resultados.append(rel)
    
    return resultados