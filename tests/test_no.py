import os
import sys

import numpy as np
import pytest

# Garante import do pacote no layout src/
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from AECPy.no import NoBase, NoGR, NoPE, NoPP, NoTE, NoTP
from AECPy.no import No_GR, No_PE, No_PP, No_TE, No_TP  # aliases legados


def test_no_base_nao_instanciavel():
    with pytest.raises(TypeError):
        NoBase([0.0, 0.0])


def test_criar_no_pp_com_coordenadas_validas():
    no = No_PP([1.0, 2.0])
    assert np.allclose(no.coor, np.array([1.0, 2.0]))
    assert no.ndim == 2
    assert no.ngdl == 3


def test_idl_e_ifl_validos():
    no = No_PP([0.0, 0.0])
    assert no.idl("ux") == 0
    assert no.ifl("fx") == 0


def test_idl_invalido_levanta_erro():
    no = No_PP([0.0, 0.0])
    with pytest.raises(ValueError):
        no.idl("uy")


def test_ifl_invalido_levanta_erro():
    no = No_PP([0.0, 0.0])
    with pytest.raises(ValueError):
        no.ifl("fy")


def test_igdl_idg_e_ifg():
    no = No_PP([0.0, 0.0])
    no.igdl = [10, 11, 12]
    assert no.idg("ux") == 10
    assert no.ifg("my") == 12


def test_deslocamentos_nulos_define_indices():
    no = No_PP([0.0, 0.0])
    no.igdl = [0, 1, 2]
    no.deslocamentos_nulos = ["uz", "ux"]
    assert set(no.deslocamentos_nulos) == {"ux", "uz"}
    assert no.idl("ux") == 0
    assert no.idl("uz") == 1
    assert no.idg("ux") == 0
    assert no.idg("uz") == 1


def test_deslocamentos_prescritos_define_indices_e_valores_em_ordem_global():
    no = No_PP([0.0, 0.0])
    no.igdl = [5, 6, 7]
    no.deslocamentos_prescritos = {"ry": 0.2, "ux": 1.0}

    assert no.deslocamentos_prescritos == {"ry": 0.2, "ux": 1.0}
    assert no.idl("ux") == 0
    assert no.idl("ry") == 2
    assert no.idg("ux") == 5
    assert no.idg("ry") == 7


def test_apoio_elastico_define_indices_e_valores_em_ordem_global():
    no = No_PP([0.0, 0.0])
    no.igdl = [3, 4, 5]
    no.apoio_elastico = {"ry": 200.0, "ux": 100.0}

    assert no.apoio_elastico == {"ry": 200.0, "ux": 100.0}
    assert no.idl("ux") == 0
    assert no.idl("ry") == 2
    assert no.idg("ux") == 3
    assert no.idg("ry") == 5


def test_check_cc_detecta_conflito_entre_nulo_e_prescrito():
    no = No_PP([0.0, 0.0])
    no.deslocamentos_nulos = ["ux"]
    with pytest.raises(Exception):
        no.deslocamentos_prescritos = {"ux": 0.0}


def test_check_cc_detecta_conflito_entre_nulo_e_elastico():
    no = No_PP([0.0, 0.0])
    no.deslocamentos_nulos = ["ux"]
    with pytest.raises(Exception):
        no.apoio_elastico = {"ux": 1e5}


def test_carga_por_dicionario_e_vetor_p():
    no = No_PP([0.0, 0.0])
    no.carga = {"fx": 10.0, "my": 2.0}
    assert bool(no.carga) is True
    assert np.allclose(no.p, np.array([10.0, 0.0, 2.0]))


def test_carga_por_lista_converte_para_dicionario():
    no = No_PP([0.0, 0.0])
    no.carga = [1.0, 2.0, 3.0]
    assert no.carga == {"fx": 1.0, "fz": 2.0, "my": 3.0}


def test_limpar_carga_com_none():
    no = No_PP([0.0, 0.0])
    no.carga = {"fx": 10.0}
    no.carga = None
    assert bool(no.carga) is False
    assert np.allclose(no.p, np.array([0.0, 0.0, 0.0]))


def test_repr_e_str_basicos():
    no = No_PP([1.0, 2.0])
    no.igdl = [0, 1, 2]
    no.carga = {"fx": 5.0}
    s = str(no)
    r = repr(no)

    assert "igdl" in s
    assert "carga" in s
    assert "NoPP([1.0, 2.0])" == r


def test_subclasses_tipos_esperados():
    assert NoPE.tipo == "PE"
    assert NoPP.tipo == "PP"
    assert NoTE.tipo == "TE"
    assert NoTP.tipo == "TP"
    assert NoGR.tipo == "GR"
