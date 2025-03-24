
import AECPy as aec


aec.No.iniciar('PP')
aec.Elemento.iniciar('PP')

# Dados de entrada [kN, cm]
L = 500.0
H = 300.0
E = 20e3

mat1 = aec.Material(E,cp=0.3,nome='aço')
sec1 = aec.SecaoRetangular(b=10.0,h=30.0,mat=mat1,nome='barra_retangular')

# Nós 
nos = [aec.No([0.0, 0.0])] * 4
print(nos)

nos[0] = aec.No([0.0, 0.0]) # type: ignore
nos[1] = aec.No([0.0, H]) # type: ignore
nos[2] = aec.No([L, H])
nos[3] = aec.No([L, 0.0])

nos[0].definir_apoio('ux','uz')
nos[3].definir_apoio('uz')
nos[1].definir_carga(fx=10)

# Elementos
els = 3*[aec.Elemento(nos[0], nos[1], sec1)]
els[0] = aec.Elemento(nos[0], nos[1], sec1)
els[1] = aec.Elemento(nos[1], nos[2], sec1)
els[2] = aec.Elemento(nos[2], nos[3], sec1)


# Análise
ngdl = len(nos) * nos[0].ngdl

ngdlF = aec.modelo.numerar_gdls(nos)
for i, no in enumerate(nos):
    print(f'{i=} : {no.igdl=}')

K, F = aec.modelo.construir_rigidez_forças(nos, els, ngdl)
U, R = aec.modelo.resolver(K, F, nos, ngdlF)

res = aec.modelo.resultados_nos(nos, U, R)

print(res )