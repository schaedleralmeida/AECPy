"""
Exemplo: pórtico plano simples
==============================
Geometria
---------
- Vão:    8 m (direção x)
- Altura: 3 m (direção z)

Seções transversais
-------------------
- Pilares: retangular 20 × 30 cm, Concreto C35
- Viga:    retangular 20 × 50 cm, Concreto C20

Carregamento
------------
- Carga distribuída na viga: 20 kN/m (vertical, para baixo)

Unidades
--------
- Comprimentos: m
- Forças:       kN

Dois modelos equivalentes são criados e comparados:
  modelo1 — definido diretamente por código
  modelo2 — carregado a partir do arquivo Excel 'portico_02.xlsx'
"""

import pandas as pd
from pathlib import Path

from AECPy import Modelo
from AECPy.secao import Material, SecaoRetangular

PASTA       = Path(__file__).parent
ARQUIVO_XLS = PASTA / "portico_02.xlsx"

# ==============================================================
# Propriedades dos materiais
# ==============================================================
# Módulos de elasticidade longitudinal (NBR 6118 — módulo inicial)
#   C20: Eci = 5600 * sqrt(20) ≈ 25 000 MPa = 25e6 kN/m²
#   C35: Eci = 5600 * sqrt(35) ≈ 33 000 MPa = 33e6 kN/m²
# Módulo de cisalhamento: G = E / (2*(1 + ν)),  ν = 0.2 → divisor 2.4
# Peso específico do concreto: 25 kN/m³

E_C20 = 25.0e6   # kN/m²
E_C35 = 33.0e6   # kN/m²
G_C20 = E_C20 / 2.4
G_C35 = E_C35 / 2.4
pe    = 25.0     # kN/m³

# ==============================================================
# Modelo 1 — definido diretamente por código
# ==============================================================

# --- Materiais ---
c35 = Material(E=E_C35, G=G_C35, pe=pe, nome="C35")
c20 = Material(E=E_C20, G=G_C20, pe=pe, nome="C20")

# --- Seções transversais (dimensões em metros) ---
sec_pilar = SecaoRetangular(c35, b=0.20, h=0.30, nome="pilar_20x30")
sec_viga  = SecaoRetangular(c20, b=0.20, h=0.50, nome="viga_20x50")

# --- Modelo ---
m1 = Modelo("PP", unidades=("m", "kN"))

# Nós: coordenadas (x, z)
m1.adicionar_no(1, [0.0, 0.0], deslocamentos_nulos=["ux", "uz", "ry"])  # base pilar esq.
m1.adicionar_no(2, [0.0, 3.0])                                           # topo pilar esq.
m1.adicionar_no(3, [8.0, 3.0])                                           # topo pilar dir.
m1.adicionar_no(4, [8.0, 0.0], deslocamentos_nulos=["ux", "uz", "ry"])  # base pilar dir.

# Elementos
m1.adicionar_elemento(1, 1, 2, sec_pilar)                           # pilar esquerdo
m1.adicionar_elemento(2, 2, 3, sec_viga, carga={"w2": -20.0})      # viga (w2 em kN/m)
m1.adicionar_elemento(3, 4, 3, sec_pilar)                           # pilar direito

m1.analisar()

# ==============================================================
# Geração do arquivo Excel  (executado apenas uma vez)
# ==============================================================

aba_materiais = pd.DataFrame([
    {"label": "C35", "E": E_C35, "G": G_C35, "pe": pe, "nome": "C35"},
    {"label": "C20", "E": E_C20, "G": G_C20, "pe": pe, "nome": "C20"},
])

aba_secoes = pd.DataFrame([
    {"label": "pilar_20x30", "tipo": "SecaoRetangular", "material": "C35", "params": "b=0.20, h=0.30"},
    {"label": "viga_20x50",  "tipo": "SecaoRetangular", "material": "C20", "params": "b=0.20, h=0.50"},
])

# Colunas 'x' e 'z' são os eixos globais do pórtico plano (NoPP.eixos_globais)
# A coluna 'deslocamentos_nulos' usa a notação de lista Python entre aspas
aba_nos = pd.DataFrame([
    {"id": 1, "x": 0.0, "z": 0.0, "deslocamentos_nulos": "['ux', 'uz', 'ry']"},
    {"id": 2, "x": 0.0, "z": 3.0},
    {"id": 3, "x": 8.0, "z": 3.0},
    {"id": 4, "x": 8.0, "z": 0.0, "deslocamentos_nulos": "['ux', 'uz', 'ry']"},
])

# A coluna 'carga' usa a notação de dicionário Python entre aspas
aba_elementos = pd.DataFrame([
    {"id": 1, "id_noI": 1, "id_noJ": 2, "sec": "pilar_20x30"},
    {"id": 2, "id_noI": 2, "id_noJ": 3, "sec": "viga_20x50",  "carga": "{'w2': -20.0}"},
    {"id": 3, "id_noI": 4, "id_noJ": 3, "sec": "pilar_20x30"},
])

with pd.ExcelWriter(ARQUIVO_XLS, engine="openpyxl") as writer:
    aba_materiais.to_excel(writer, sheet_name="materiais", index=False)
    aba_secoes.to_excel(writer,    sheet_name="secoes",    index=False)
    aba_nos.to_excel(writer,       sheet_name="nos",       index=False)
    aba_elementos.to_excel(writer, sheet_name="elementos", index=False)

print(f"Arquivo Excel gerado: {ARQUIVO_XLS.name}\n")

# ==============================================================
# Modelo 2 — carregado a partir do arquivo Excel
# ==============================================================

m2 = Modelo("PP", unidades=("m", "kN"))
m2.carregar_excel(ARQUIVO_XLS)
m2.analisar()

# ==============================================================
# Resultados
# ==============================================================

print("=== Modelo 1 (código) ===")
print(m1)
print()
print(m1.tabela_nos())

print("\n=== Modelo 2 (Excel) ===")
print(m2)
print()
print(m2.tabela_nos())
