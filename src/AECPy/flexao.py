'''
Rotinas relativas ao comportamento à flexão do elmenento
em torno dos eixos princiapais centrais 2 e 3 da seção transversal
'''

import numpy as np

#-----------------------------------------------------------
#          flexão em torno do eixo local 3
#        -- vigas, grelhas, pórticos planos e espaciais --
#-----------------------------------------------------------

def calc_Kb3(S,L):
    '''
    Cálcul a matriz de rigidez local para a flexão em torno do eixo local 3
    '''
    r1 = 12*S.EI3/L**3
    r2 =  6*S.EI3/L**2
    r3 =  4*S.EI3/L

    return np.array([[ r1 , r2  ,-r1, r2   ],
                     [ r2 , r3  ,-r2, r3/2 ],
                     [-r1 ,-r2  , r1,-r2   ],
                     [ r2 , r3/2,-r2, r3   ]])

def rel_d_b3(xi,S,L,db3):
    '''Deslocamentos e esforços pro flexão em torno do eixo 3
    calculados em função dos deslocamentos nodais'''
    
    'deslocamentos nodais'
    u2i,r3i,u2j,r3j = db3
    
    '''Coordenada adimensional xi '''
    xi2 = xi**2
    xi3 = xi**3

    'Rigidez da seção transversal'
    EI = S.EI3
    
    u_2 = ((2*xi3 - 3*xi2 + 1 ) *u2i
           + (xi3 - 2*xi2 + xi)   *r3i*L
           + (-2*xi3 + 3*xi2)     *u2j
           + (xi3 -xi2)           *r3j*L
           )
    
    rot_3 = ((6*xi2 - 6*xi )    *u2i/L
             + (3*xi2 -4*xi +1)   *r3i
             + (-6*xi2 + 6*xi)    *u2j/L
             + (3*xi2 - 2*xi)     *r3j
             )
    
    V2 = 12*EI/L**3 * ( (u2j - u2i) - L*(r3j+r3i)/2 )  + 0*xi
    
    M3 = (6*EI)/(L**2) * (2*xi-1)*(u2i-u2j) +(2*EI)/L * ( (3*xi-2)*r3i +(3*xi-1)*r3j)
    
    return {'u2':u_2, 'r3':rot_3, 'V2':V2, 'M3':M3}


def rep_w2(L,w2):
    '''Reação de engastamento perfeito da força distribuída
     no sentido do eixo local 2 - associado à flexão em torno do eixo 3'''
     
    re_b3 = np.zeros(4)
    Dw2 = w2[1] - w2[0]
    re_b3[0] = -( w2[0]*L/2 + 3*Dw2*L/20 )
    re_b3[1] = -( w2[0]*L**2/12 + Dw2*L**2/30)
    re_b3[2] = -( w2[0]*L/2 + 7*Dw2*L/20 )
    re_b3[3] =  ( w2[0]*L**2/12 + Dw2*L**2/20)
    return re_b3


def rel_w2(xi,S,L,w2):
    '''Deslocamentos e esforços no ponto de coordenada adimensional xi
     no elemento calculados em função dos efeitos locais da força
     distribuída no sentido do eixo local 2
     - associado à flexão em torno do eixo 3 (no plano 1-2)'''
     
    xi2 = xi**2
    xi3 = xi**3
    xi4 = xi**4
    xi5 = xi**5
    
    'Rigidez da seção transversal'
    EI = S.EI3
    
    #carregamento com variação linear em todo o elemento
    Dw2 = w2[1]-w2[0]
    u_2 =  ( (w2[0]*L**4)/(24*EI) * (xi4 - 2*xi3 + xi2 )
             + (Dw2*L**4)/(120*EI) * (xi5 - 3*xi3 + 2*xi2 ) )
    rot_3 = ( (w2[0]*L**3)/(12*EI) * (2*xi3 - 3*xi2 + xi)
              + (Dw2*L**3)/(120*EI) * (5*xi4 - 9*xi2 + 4*xi) )
    V2 = -w2[0]*L/2 *(2*xi-1) + Dw2*L/20*(3-10*xi2)
    M3 = w2[0]*L**2/12 *(6*xi2 - 6*xi +1) + Dw2*L**2/60*(10*xi3-9*xi+2)
    
    return {'u2':u_2, 'r3':rot_3, 'V2':V2, 'M3':M3}


def rep_T2(L, S, dT):
    '''Reações de engastamento perfeito pela variação de temperatura na direção do eixo 2 da seção transversal da barra'''
    #obs: dT é a variação da temperatura do centroide à face superior ( perpendicular ao eixo 2) 
    return dT * S.mat.cdt * S.mat.E * S.W3 * np.array([0,-1,0,1])
    
def rel_T2(xi,S,L,dT):
    '''Deslocamentos e esforços no ponto de coordenada adimensional xi
     no elemento calculados em função dos efeitos locais da variação de temperatura na direção do eixo 2 da seção transversal da barra'''
    #obs: dT é a variação da temperatura do centroide à face superior ( perpendicular ao eixo 2)  
    return {'M3': dT * S.mat.cdt * S.mat.E * S.W3 * xi**0}

#-----------------------------------------------------------
#          flexão em torno do eixo local 2
#             -- pórticos espaciais -- 
#-----------------------------------------------------------

def calc_Kb2(S,L):
    '''
    Cálcul a matriz de rigidez local para a flexão em torno do eixo local 2
    '''        

    r1 = 12*S.EI2/L**3
    r2 =  6*S.EI2/L**2
    r3 =  4*S.EI2/L
        
    return np.array([[ r1 ,-r2  ,-r1,-r2   ],
                     [-r2 , r3  , r2, r3/2 ],
                     [-r1 , r2  , r1, r2   ],
                     [-r2 , r3/2, r2, r3   ]])


def rel_d_b2(xi,S,L,db2):
    '''Deslocamentos e esforços pro flexão em torno do eixo 2
    calculados em função dos deslocamentos nodais'''
    
    'deslocamentos nodais'
    u3i,r2i,u3j,r2j = db2
    
    '''Coordenada adimensional xi '''
    xi2 = xi**2
    xi3 = xi**3

    'Rigidez da seção transversal'
    EI = S.EI2
    
    u_3 = (  (2*xi3 - 3*xi2 + 1 ) *u3i
           - (xi3 - 2*xi2 + xi)   *r2i*L
           + (-2*xi3 + 3*xi2)     *u3j
           - (xi3 -xi2)           *r2j*L
           )
    
    rot_2 = (  (6*xi2 - 6*xi )    *u3i/L
             - (3*xi2 -4*xi +1)   *r2i
             + (-6*xi2 + 6*xi)    *u3j/L
             - (3*xi2 - 2*xi)     *r2j
             )
    
    V3 = 12*EI/L**3 * ( (u3j - u3i) + L*(r2j+r2i)/2 ) + 0*xi
    
    M2 = (6*EI)/(L**2) * (2*xi-1)*(u3i-u3j) -(2*EI)/L * ( (3*xi-2)*r2i +(3*xi-1)*r2j)
    
    return {'u3':u_3, 'r2':rot_2, 'V3':V3, 'M2':M2}

def rep_w3(L,w3):
    '''Reação de engastamento perfeito da força distribuída
     no sentido do eixo local 3 - associado à flexão em torno do eixo 2'''
    
    re_b2 = np.zeros(4)
    Dw3 = w3[1] - w3[0]
    re_b2[0] = -( w3[0]*L/2 + 3*Dw3*L/20 )
    re_b2[1] =  ( w3[0]*L**2/12 + Dw3*L**2/30)
    re_b2[2] = -( w3[0]*L/2 + 7*Dw3*L/20 )
    re_b2[3] = -( w3[0]*L**2/12 + Dw3*L**2/20)
    return re_b2


def rel_w3(xi,S,L,w3):
    '''Deslocamentos e esforços no ponto de coordenada adimensional xi
     no elemento calculados em função dos efeitos locais da força
     distribuída no sentido do eixo local 3
     - associado à flexão em torno do eixo 2 (no plano 1-3)'''
     
    xi2 = xi**2
    xi3 = xi**3
    xi4 = xi**4
    xi5 = xi**5
    
    'Rigidez da seção transversal'
    EI = S.EI2
    
    #carregamento com variação linear em todo o elemento
    Dw3 = w3[1]-w3[0]
    u_3 =  ( (w3[0]*L**4)/(24*EI) * (xi4 - 2*xi3 + xi2 )
             + (Dw3*L**4)/(120*EI) * (xi5 - 3*xi3 + 2*xi2 ) )
    rot_2 = - ( (w3[0]*L**3)/(12*EI) * (2*xi3 - 3*xi2 + xi)
              - (Dw3*L**3)/(120*EI) * (5*xi4 - 9*xi2 + 4*xi) )
    V3 = -w3[0]*L/2 *(2*xi-1) + Dw3*L/20*(3-10*xi2)
    M2 = w3[0]*L**2/12 *(6*xi2 - 6*xi +1) + Dw3*L**2/60*(10*xi3-9*xi+2)
    
    return {'u3':u_3, 'r2':rot_2, 'V3':V3, 'M2':M2}


def rep_T3(L, S, dT):
    '''Reações de engastamento perfeito pela variação de temperatura na direção do eixo 3 da seção transversal da barra'''
    #obs: dT é a variação da temperatura do centroide à face superior ( perpendicular ao eixo 3) 
    return dT * S.mat.cdt * S.mat.E * S.W2 * np.array([0,1,0,-1])
    
def rel_T3(xi,S,L,dT):
    '''Deslocamentos e esforços no ponto de coordenada adimensional xi
     no elemento calculados em função dos efeitos locais da variação de temperatura na direção do eixo 3 da seção transversal da barra'''
     #obs: dT é a variação da temperatura do centroide à face superior ( perpendicular ao eixo 3) 
    return {'M2': dT * S.mat.cdt * S.mat.E * S.W2 * xi**0}