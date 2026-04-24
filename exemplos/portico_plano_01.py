'''
Exemplo 4.13 do livro Matrix Structural Analysis 2ed
'''
#%%
import numpy as np
from math import sin, cos, radians

import AECPy as aec
from AECPy.no import NoPP as No
from AECPy.elemento import ElementoPP as Elemento
from AECPy.secao import Secao, Material

#unidades: [kN, mm]

#%%
#materiais:
mat1 = Material(200, nome='meu_material')

#seções transversais (adotando I2=J=1 em todas):
sec_ab  = Secao(mat1, 6.e3, I3=200.e6, nome='Sec_AB')
sec_bc  = Secao(mat1, 4.e3, I3= 50.e6, nome='Sec_BC')

#%%
# Metodologia 1: análise procedural (funções de módulo)
nos = []
nos.append(No([0.0, 5.0e3]))
nos.append(No([8.0e3, 5.0e3]))
nos.append(No([8.0e3, 0.0]))
nos[0].deslocamentos_nulos = list(No.gdls_globais)
nos[2].deslocamentos_nulos = list(No.gdls_globais)

ang = radians(-45)
nos[1].carga = {'fx': 100 * cos(ang), 'fz': 100 * sin(ang), 'my': -50.e3}

els = []
els.append(Elemento(nos[0], nos[1], sec_ab))
els.append(Elemento(nos[1], nos[2], sec_bc))

ngdl = len(nos) * nos[0].ngdl
ngdlF = aec.modelo.numerar_gdls(nos)
K, F = aec.modelo.construir_SEL(nos, els, ngdl)
U_proc, R_proc = aec.modelo.resolver_SEL(K, F, nos, ngdlF)
res_nos_proc_list = aec.modelo.resultados_nos(nos, U_proc, R_proc)
res_nos_proc = {i: s for i, s in enumerate(res_nos_proc_list)}
res_els_proc = aec.modelo.resultados_elementos(els, U_proc)

print('\n=== Metodologia procedural ===')
print('Deslocamentos e Reações na direção x para os nós 1 e 2:')
print(np.array([[res_nos_proc[1]['ux'], res_nos_proc[1]['Rfx']],
				[res_nos_proc[2]['ux'], res_nos_proc[2]['Rfx']]]))
print(f"M3 no centro do elemento 0: {res_els_proc[0]['M3'][2]:.2f} kN.mm")

#%%
# Metodologia 2: análise com objeto Modelo
modelo = aec.modelo.Modelo('PP', unidades=('mm', 'kN'))

modelo.adicionar_no(1, [0.0, 5.0e3], deslocamentos_nulos=['ux', 'uz', 'ry'])
modelo.adicionar_no(2, [8.0e3, 5.0e3], carga={'fx': 100 * cos(ang), 'fz': 100 * sin(ang), 'my': -50.e3})
modelo.adicionar_no(3, [8.0e3, 0.0], deslocamentos_nulos=['ux', 'uz', 'ry'])

modelo.adicionar_elemento(1, 1, 2, sec_ab)
modelo.adicionar_elemento(2, 2, 3, sec_bc)

#%%
modelo.analisar(npts=5)
res_nos_obj = modelo.res_nos
res_els_obj = modelo.res_elementos
assert modelo.U is not None and modelo.R is not None
assert res_nos_obj is not None and res_els_obj is not None

print('\n=== Metodologia com objeto Modelo ===')
print('Deslocamentos e Reações na direção x para os nós 2 e 3:')
print(np.array([[res_nos_obj[2]['ux'], res_nos_obj[2]['Rfx']],
				[res_nos_obj[3]['ux'], res_nos_obj[3]['Rfx']]]))
print(f"M3 no centro do elemento 1: {res_els_obj[1]['M3'][2]:.2f} kN.mm")

#%%
# Comparação entre metodologias
ok_u = np.allclose(U_proc, modelo.U)
ok_r = np.allclose(R_proc, modelo.R)
ok_m3 = np.isclose(res_els_proc[0]['M3'][2], res_els_obj[1]['M3'][2])

print('\n=== Comparação ===')
print(f'Deslocamentos iguais? {ok_u}')
print(f'Reações iguais? {ok_r}')
print(f'Momento M3 no centro do elemento correspondente igual? {ok_m3}')

print('\nTabela de resultados (objeto Modelo):')
print(modelo.tabela_nos())
print('\nTabela do elemento 1 (objeto Modelo):')
print(modelo.tabela_elemento(1, ['M3', 'V2']))

fig = modelo.diagramas_elemento(1, ['u2', 'M3', 'V2'])
fig_modelo = modelo.visualizar()

#%%
# Escrita dos resultados em planilha
aec.planilhas.salvar_excel('portico_plano_01.xlsx', modelo)
print('\nResultados salvos em portico_plano_01.xlsx')

#%%
# Reconstrução do modelo a partir da planilha gerada
modelo_reconstruido = aec.modelo.Modelo('PP', unidades=('mm', 'kN'))
aec.planilhas.carregar_excel('portico_plano_01.xlsx', modelo_reconstruido)
modelo_reconstruido.analisar(npts=5)
print('\n=== Modelo reconstruído a partir da planilha ===')
print(modelo_reconstruido)
print(modelo_reconstruido.tabela_nos())

#%%
# Verificação: resultados do modelo reconstruído são iguais ao original?
assert modelo.U is not None and modelo_reconstruido.U is not None
assert modelo.R is not None and modelo_reconstruido.R is not None
ok_u   = np.allclose(modelo.U, modelo_reconstruido.U)
ok_r   = np.allclose(modelo.R, modelo_reconstruido.R)

res_els_rec = modelo_reconstruido.res_elementos
assert res_els_rec is not None
ok_els = all(
    np.allclose(res_els_obj[id_el][g], res_els_rec[id_el][g])
    for id_el in res_els_obj
    for g in res_els_obj[id_el]
)

print('\n=== Verificação modelo original × reconstruído ===')
print(f'Deslocamentos iguais?        {ok_u}')
print(f'Reações iguais?              {ok_r}')
print(f'Esforços nos elementos iguais? {ok_els}')