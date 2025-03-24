"""
Rotinas relativas ao comportamento `axial` do elmenento:
compressão/tração e torção
"""

import numpy as np

# -----------------------------------------------------------
#                    axial
# -----------------------------------------------------------


def calc_Ka(S, L):
    """
    Calcula a matriz de rigidez local para a deformação axial
    """

    ra = S.EA / L
    return np.array([[ra, -ra], [-ra, ra]])


def rel_d_a(xi, S, L, da):
    """Deslocamentos e esforços por deformação axial
    calculados em função dos deslocamentos nodais"""

    u1i, u1j = da
    u_1 = u1i * (1 - xi) + u1j * xi
    N = S.EA / L * (u1j - u1i) + 0 * xi
    return {"u1": u_1, "N": N}


def rep_w1(L, w1):
    """Reação de engastamento perfeito da força distribuída
    no sentido do eixo local 1 - associado à deformação axial"""

    re_a = np.zeros(2)
    Dw1 = w1[1] - w1[0]
    re_a[0] = -w1[0] * L / 2 - Dw1 * L / 6
    re_a[1] = -w1[0] * L / 2 - Dw1 * L / 3
    return re_a


def rel_w1(xi, S, L, w1):
    """Deslocamentos e esforços no ponto de coordenada adimensional xi
    no elemento calculados em função dos efeitos locais da força
    distribuída no sentido do eixo local 1 - associado à deformação axial"""

    xi2 = xi**2
    xi3 = xi**3

    "Rigidez da seção transversal"
    EA = S.EA

    Dw1 = w1[1] - w1[0]
    u_1 = (w1[0] * L**2) / (2 * EA) * (xi - xi2) + (Dw1 * L**2 / (6 * EA)) * (
        xi - xi3
    )
    N = (w1[0] * L) / (2) * (1 - 2 * xi) + (Dw1 * L / 6) * (1 - 3 * xi2)

    return {"u1": u_1, "N": N}


def rep_T(L, S, dT):
    """Reações de engastamento perfeito pela variação uniforme de temperatura na barra"""
    return dT * S.mat.cdt * S.EA * np.array([1, -1])


def rel_T(xi, S, L, dT):
    """Deslocamentos e esforços no ponto de coordenada adimensional xi
    no elemento calculados em função dos efeitos locais da variação de temperatura"""
    return {"N": -S.mat.cdt * dT * S.EA * xi**0}


# -----------------------------------------------------------
#                    torção
# -----------------------------------------------------------
def calc_Kt(S, L):
    """
    Cálcul a matriz de rigidez local para a torção pura
    """

    rt = S.GJ / L
    return np.array([[rt, -rt], [-rt, rt]])


def rel_d_t(xi, S, L, dt):
    """Deslocamentos e esforços por torção
    calculados em função dos deslocamentos nodais"""

    r1i, r1j = dt
    rot_1 = r1i * (1 - xi) + r1j * xi
    Mt = S.GJ / L * (r1j - r1i) + 0 * xi
    return {"r1": rot_1, "Mt": Mt}
