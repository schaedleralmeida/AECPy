'''
Exemplo 4.13 do livro Matrix Structural Analysis 2ed
'''
import numpy as np
from math import sin, cos, radians

import AECPy as aec
from AECPy.no import NoPP as No
from AECPy.elemento import ElementoPP as Elemento
from AECPy.secao import Secao, Material

#unidades: [kN, mm]

#materiais:
mat1 = Material(200, nome='meu_material')

#seções transversais (adotando I2=J=1 em todas):
sec_ab  = Secao(mat1, 6.e3, I3=200.e6, nome='Sec_AB')
sec_bc  = Secao(mat1, 4.e3, I3= 50.e6, nome='Sec_BC')

#nós:
nos = []
nos.append(No([0.0, 5.0e3]))
nos.append(No([8.0e3, 5.0e3]))
nos.append(No([8.0e3, 0.0]))
nos[0].deslocamentos_nulos = list(No.gdls_globais)
nos[2].deslocamentos_nulos = list(No.gdls_globais)
ang = radians(-45)
nos[1].carga = {'fx': 100 * cos(ang), 'fz': 100 * sin(ang), 'my': -50.e3}

#elementos:
els = []
els.append(Elemento(nos[0], nos[1], sec_ab))
els.append(Elemento(nos[1], nos[2], sec_bc))


# Análise
ngdl = len(nos) * nos[0].ngdl

ngdlF = aec.modelo.numerar_gdls(nos)
# for i, no in enumerate(nos):
#     print(f'{i=} : {no.igdl=}')

K, F = aec.modelo.construir_SEL(nos, els, ngdl)
# with np.printoptions(precision=0,linewidth=100,suppress=True):
#     print(K)
#     print(F)
U, R = aec.modelo.resolver_SEL(K, F, nos, ngdlF)

res_nos = aec.modelo.resultados_nos(nos, U, R)

print('\n Deslocamentos e Reações na direção x par os nós 1 e 2:')
print( res_nos.loc[[1,2],['ux','Rfx']])


res_els = aec.modelo.resultados_elementos(els, U)
print(f'\n O momento M3 no centro do elemento 0 é {res_els[0]["M3"][2]:.2f} kNmm')

print('\n O momento e força cortante no elemento 0 são:')
uni = aec.unidades.Conversor('mm', 'kN')
tab = aec.graficos.tabela_rel(res_els[0], ['M3', 'V2'], uni)
print(tab)


print('\n O diagrama de deslocamento transversal, momento e força cortante no elemento 0 são:')
fig = aec.graficos.diagramas(res_els[0], ['u2', 'M3', 'V2'], uni, eltag=0)



print('Representação do modelo:')
fig_modelo = aec.graficos.modelo_2d(nos, els)