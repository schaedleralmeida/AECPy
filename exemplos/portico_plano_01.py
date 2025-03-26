'''
Exemplo 4.13 do livro Matrix Structural Analysis 2ed
'''
import numpy as np
from math import sin, cos, radians

import sys
import os
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.realpath(f"{dir_path}/../src/"))

import AECPy as aec

aec.No.iniciar('PP')
aec.Elemento.iniciar('PP')

#unidades: [kN, mm]

#materiais:
mat1 = aec.Material(200,0.0,nome='meu_material')

#seções transversais (adotando I2=J=1 em todas):
sec_ab  = aec.Secao(mat1, 6.e3, I3=200.e6, nome='Sec_AB') 
sec_bc  = aec.Secao(mat1, 4.e3, I3= 50.e6, nome='Sec_BC')

#nós:
nos = [aec.No([0.0,0.0])] * 3
nos[0] = aec.No([0.0, 5.0e3])
nos[1] = aec.No([8.0e3, 5.0e3])
nos[2] = aec.No([8.0e3, 0.0])
nos[0].definir_apoio('todos')
nos[2].definir_apoio('todos')
ang = radians(-45)
nos[1].definir_carga(fx=( 100 * cos(ang)) , fz = ( 100 * sin(ang)), my = -50.e3)

#elementos:
els = [aec.Elemento(nos[0], nos[1], sec_ab)] * 2
els[0] = aec.Elemento(nos[0], nos[1], sec_ab)
els[1] = aec.Elemento(nos[1], nos[2], sec_bc)


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





#exeplo de resolução usado a classe Modelo
print( 3*'\n', 15*'=', 'Exemplo com a classe Modelo', 15*'=', 3*'\n')
uni = aec.unidades.Conversor('m', 'kN')

#materiais:
mat1 = aec.Material(200*uni.em('GPa'),0.0,nome='meu_material')

#seções transversais (adotando I2=J=1 em todas):
sec_ab  = aec.Secao(mat1, 6.e3*uni.em('mm2'), I3=200.e6*uni.em('mm4'), nome='Sec_AB') 
sec_bc  = aec.Secao(mat1, 4.e3*uni.em('mm2'), I3= 50.e6*uni.em('mm4'), nome='Sec_BC')


#nós:
nos = [aec.No([0.0,0.0])] * 3
nos[0] = aec.No([0.0, 5.0])
nos[1] = aec.No([8.0, 5.0])
nos[2] = aec.No([8.0, 0.0])
nos[0].definir_apoio('todos')
nos[2].definir_apoio('todos')
ang = radians(-45)
nos[1].definir_carga(fx=( 100 * cos(ang)) , fz = ( 100 * sin(ang)), my = -50.e3)

#elementos:
els = [aec.Elemento(nos[0], nos[1], sec_ab)] * 2
els[0] = aec.Elemento(nos[0], nos[1], sec_ab)
els[1] = aec.Elemento(nos[1], nos[2], sec_bc)
estr = aec.Modelo(nos, els,uni)
estr.npts_res_els = 10
estr.resolver()
diagrama_el0 = estr.diagramas(0, ['N', 'M3', 'V2'])

print( estr.res_nos.loc[[1,2],['ux','Rfx']])
fig_diagrama  = aec.graficos.diagramas( estr.res_els[0], ['u2', 'M3', 'V2'])
fig_modelo_2 = estr.figura()
fig_modelo_2.show()
