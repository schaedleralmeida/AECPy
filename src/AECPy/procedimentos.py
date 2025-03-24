"""
Contém procedimentos básicos de métodos matriciais que serão usados
para análise estrutural pelo Método da Resistência Direta (MRD) no AECPy
"""

from math import cos, radians, sin, sqrt

import numpy as np


def calc_L_u(pI, pJ):
    """
    Calcula o comprimento e o vetor unitário que dá a direção de um segmento de reta

    Parameters
    ----------
    pI,pJ: numpy.ndarray
        Coordenadas da posição inicial e final do segmento de reta
        pI.shape = pJ.shape = (n,), onde `n` é o número de coordenadas que define a posição

    Returns
    -------
    L: float
        Comprimento do segmento de reta
    u: numpy.ndarray
        Vetor unitário apontando de pI para pJ. u.shape = (n,)
    """
    rij = pJ - pI
    L = np.linalg.norm(rij)
    u = rij / L
    return L, u


def reunir(Ag, ig) -> np.ndarray:
    """
    Gather: Recolhe quantidade referente aos gdl do nó/elemento de um array global (Ag)
    e armazena em um array local (Al) do elemento.

    Parameteres
        -----------
        Ag: numpy.ndarray
            Array global com os valores a serem obtidos. (Ag.ndim = 1 ou 2)
        ig: iterable (list, numpy.ndarray) contendo int
            Índices globais referentes aos gdl de interesse para o nó/elemento
        Returns
        -------
        Al: numpy.ndarray
            Array com quantidades referentes ao nó/elemento
            Se Ag.ndim = 1, Al.shape = (len(ig),)
            Se Ag.ndim = 2, Al.shape = (len(ig),len(ig))
    """
    if not isinstance(ig, (list, np.ndarray)):
        ig = list(ig)

    # caso em que Ag é um array unidimensional (vetor)
    if Ag.ndim == 1:
        Al = np.array(Ag[ig])

    # caso em que Ag é um array bidimensional (matriz)
    elif Ag.ndim == 2:
        nig = len(ig)  # número índices globais
        Al = np.ndarray(
            (nig, nig)
        )  # array local com as informações recolhidas

        for il in range(nig):  # preenchendo cada linha de Al
            Al[il, :] = Ag[ig[il], ig]
    else:
        raise ValueError("Ag.ndim deve ser 1 ou 2")

    return Al


def espalhar(Al, Ag, ig) -> None:
    """
    Scatter: Soma quantidade referente aos gdl do nó/elemento, armazenadas em um array local (Al),
    a um array global (Ag)

    Parameteres
    -----------
    Al: numpy.ndarray
        Array local com os valores a serem espalhados. (Al.ndim = 1 ou 2)
    Ag: numpy.ndarray
        Array global onde seão somados valores a serem espalhados. (Ag.ndim = Al.ndim)
    ig: iterable (list, tuple, numpy.ndarray) contendo int
        Índices globais referentes aos gdl de interesse para o nó/elemento
    """
    if not isinstance(ig, (list, np.ndarray)):
        ig = list(ig)

    # caso em que Ag e Al são arrays unidimensionais (vetores)
    if Ag.ndim == 1 and Al.ndim == 1:
        Ag[ig] += Al[:]

    # caso em que Ag e Al são arrays bidimensionais (matrizes)
    elif Ag.ndim == 2 and Al.ndim == 2:
        nig = len(ig)  # número índices globais

        for il in range(nig):  # preenchendo cada linha de Al
            Ag[ig[il], ig] += Al[il, :]

    return


def transf_coord(A, R, inv=False) -> np.ndarray:
    """Operação de transformação de quantidades entre sistemas de coordenadas

    Parameteres
    -----------
    A: numpy.ndarray
        Array com as quantidades para transformação de coordenadas.
        A.shape = (nA,) ou (nA,nA), respeitando nA % nR = 0
    R: numpy.ndarray
        Matriz ortogonal que define a transformação de coordenadas.
        R.shape = (nR,nR)
    inv: bool
        Indica se deve ser realizada a operação inversa à definida por `R`

    Returns
    -------
    B: numpy.ndarray
        Array com as quantidades após a transformação de coordenadas.
        B.shape = (nA,) ou (nA,nA)
    """

    if A.shape[0] % R.shape[0] != 0:
        raise ValueError("A transformação de `A` por `R` não é compatível")
    if R.shape[0] != R.shape[1]:
        raise ValueError("`A` deve ser uma matriz quadrada")
    if A.ndim > 2 or (A.ndim == 2 and A.shape[0] != A.shape[1]):
        raise ValueError(
            "`A` deve ser um array unidmensional ou bidimensional quadrada"
        )

    s = R.shape[0]  # dimensão de R
    n = (
        A.shape[0] // s
    )  # Número de "blocos diagonal com a matriz R para a rotação"
    B = np.zeros(A.shape)  # Resultado da transformação

    # matriz usada na transformação
    # (R ou transposta de R se indicada transformação inversa)
    RR = np.transpose(R) if inv else R

    # produto diag(R) A
    for i in range(n):
        i1 = i * s
        i2 = i1 + s
        B[i1:i2] = RR @ A[i1:i2]

    # para o caso em que `A` é uma matriz
    # produto diag(R) A diag(Rt)
    if A.ndim == 2:
        Rt = R if inv else np.transpose(RR)  # Transposta de R

        for i in range(n):
            i1 = i * s
            i2 = i1 + s
            B[:, i1:i2] = B[:, i1:i2] @ Rt

    return B


def R3D(e1) -> np.ndarray:
    """
    Calcula amatriz de transformação entre as coordenadas globais x-y-z
    e as coordenadas locais 1-2-3 em um problema tridimensional (3D)
    Para um vetor v = [vx, vy, vz], temos a transformação para os eixos locais
    dada por [v1, v2, v3]^T = [R] [vx, vy, vz ]^T

    Parameters
    ----------
    e1: numpy.ndarray
        Vetor unitário na direção do início ao fim do elemento (segmento de reta)
        e1.shape = (3,)

    Returns
    -------
    R : numpy.ndarray
        Matriz ortogonal de transformação de coordenadas Globais para locais;
        R.shape = (3,3)
        R=[[e1],[e2],[e3]], onde ei é o vetor unitário na direção do eixo local i=1,2,3
    """

    if e1.shape != (3,):
        raise ValueError("e1 incorreto")

    tol_verticalidade = 1.0e-3
    vertical = sqrt(1 - e1[2] ** 2) < tol_verticalidade

    if vertical:
        e2 = np.array([1, 0, 0])
        e3 = np.cross(e1, e2)

    else:  # elemento não é vertical
        aux3 = np.cross(e1, np.array([0, 0, 1]))
        e3 = aux3 / np.linalg.norm(aux3)
        e2 = np.cross(e3, e1)

    R = np.array([e1, e2, e3])

    return R


def R3D_mod_ang(R, ang):
    """
    Calcula a matriz de transformação global - local
    modificada pela rotação dos eixos locais em torno de e1 pelo ângulo `ang`

    Parameters
    ----------
    R: numpy.ndarray
        Matriz de transformação global-local padrão (R.shape = (3,3))
    ang: float
        Ângulo de rotação dos eixos locais 2-3 em torno do eixo 1
        em graus e com sentido positivo dado pela regra da mão direita

    Retorna
    -------
    Rmod: numpy.ndarray
        Matriz de transformação global-local com eixos 2-3 girados em
        relação à direção padrão (Rang.shape = (3,3))
    """
    if ang == 0:
        return R
    "rotação dos eixos locais 2-3"
    ang_rad = radians(ang)
    s = sin(ang_rad)
    c = cos(ang_rad)
    Rang = np.array([[1, 0, 0], [0, c, s], [0, -s, c]])
    Rmod = Rang @ R
    return Rmod


def R3D_u(e1, u):
    """
    Calcula amatriz de transformação entre as coordenadas globais x-y-z
    e as coordenadas locais 1-2-3 em um problema tridimensional (3D)
    considerando que o eixo local 2 está no plano formado por `e1 e `u`
    apontando na direção positiva de `u`

    Para um vetor v = [vx, vy, vz], temos a transformação para os eixos locais
    dada por [v1, v2, v3]^T = [R] [vx, vy, vz ]^T

    Parameters
    ----------
    e1: numpy.ndarray
        Vetor unitário na direção do início ao fim do elemento (segmento de reta)
        e1.shape = (3,)
    u: numpy.ndarray
        Vetor unitário usado para definir direção e sentido de e2
        e2 está no plano e1-u e tem projeção positiva em u (e2@u>0)
        u.shape = (3,)

    Returns
    -------
    R : numpy.ndarray
        Matriz ortogonal de transformação de coordenadas Globais para locais;
        R.shape = (3,3)
        R=[[e1],[e2],[e3]], onde ei é o vetor unitário na direção do eixo local i=1,2,3
    """
    if e1.shape != (3,):
        raise ValueError("e1 incorreto")
    if u.shape != (3,):
        raise ValueError("u incorreto")

    tol_alinhamento = 1.0e-3
    if 1 - e1 @ u < tol_alinhamento:
        raise ValueError("e1 e u estão alinhados")

    aux3 = np.cross(e1, u)
    e3 = aux3 / np.linalg.norm(aux3)
    e2 = np.cross(e3, e1)

    R = np.array([e1, e2, e3])

    return R


def reordenar_array(A, ord_lin, ord_col=None):
    """
    Reordena as linhas (e colunas) de uma array

    Parametres
    ---------
    A: numpy.ndarray
        array que será reordenado. A.ndim = 1 ou 2
    ord_lin: iterable
        ordem que as linhas de A devem aparacer após o reordenamento
    ord_col: iterable (optional)
        ordem em que as colunas de A devem aparecer no array reordenado
    """
    A = A[ord_lin]
    if A.ndim == 2:
        A = A[:, ord_col]
    elif A.ndim > 2:
        raise ValueError("A.dim > 2 não é suportado")
    return A


def check_xi(xi):
    """Verificação do valor de xi - coordenada adimensional entre 0 e 1"""
    if isinstance(xi, np.ndarray):
        err_xi = np.any(xi < 0) or np.any(xi > 1)
    else:
        err_xi = xi < 0 or xi > 1

    if err_xi:
        raise ValueError("xi fora do intervalo [0,1]")

    return


def igdl_FS(ngdl, ilr):
    """
    Define a numeração do índice global dos gdl dos nós,
    separando os gdl livres (F) e os gdl restritos/"de apoio" (S)

    Parameters
    ----------
    ngdl: int
        número de gdl por nó
    ilr: list
        lista com os índices locais dos gdl restritos de cada nó
        exemplo: [ [], [0,2], [], [1]] indica um modelo com 4 nós, onde
        os deslocamentos são nulos nos gdl de índices locais 0 e 2 do 2º
        e índice 1 do 4º nó.

    Returns
    -------
    igdl: numpy.ndarray
        índices globais dos gdls dos nós.
        igdl[n,l] é o índice global do gdl local `l` do nó `n`
    (nF,nS): tuple
        nF é o número de gdl livres e nS é o número de gdl restritos (apoios)
    """

    # número total de nós do modelo
    nnos = len(ilr)
    # número total de gdl do modelo
    ngdl_total = nnos * ngdl
    # número total de graus de liberdade 'restritos'/'de apoio'
    nS = 0
    for ii in ilr:
        nS += len(ii)
        # if not all([isinstance(i,int) and 0>=i>=ngdl for i in ii]) or len(ii)>:
        #     raise ValueError('ilr com valor incorreto')
    # número total de graus de liberdade 'livres'
    nF = ngdl_total - nS

    # array auxiliar como o índice local dos gdl do nó
    il = np.array([i for i in range(ngdl)], dtype=int)

    # array com o índice global dos gdl dos nós
    igdl = np.zeros((nnos, ngdl), dtype=int)

    # definção dos índices globais dos gdl dos nós
    inc_F = 0
    inc_S = nF
    for no in range(nnos):
        n = len(ilr[no])  # número de gdl restritos no nó
        if n == 0:  # nós em gdl restritos
            igdl[no, :] = il + inc_F
            inc_F += ngdl
        else:  # nós com `n` gdl restritos
            for i in range(ngdl):  # índice local do gdl
                if i in ilr[no]:  # está entre os restritos no nó
                    igdl[no, i] = inc_S
                    inc_S += 1
                else:  # não está entre os restritos no nó
                    igdl[no, i] = inc_F
                    inc_F += 1
    return igdl, (nF, nS)


# -------------------------------------------------------------------------------
#  sem uso
# -------------------------------------------------------------------------------
def expandir_array(A, index):
    m = A.shape[0]
    n = len(index) + m

    jj = [i not in index for i in range(n)]

    ii = []
    for j in range(n):
        if jj[j]:
            ii.append(j)
    print(ii)

    if A.ndim == 1:
        eA = np.zeros(n)
        print(eA.shape, A.shape, ii)
        eA[ii] = A[:]

    elif A.ndim == 2 and A.shape[0] == A.shape[1]:
        eA = np.zeros((n, n))
        for j, i in enumerate(ii):
            eA[i, ii] = A[j, :]
    else:
        raise ValueError("o formato de A não é suportado")

    return eA
