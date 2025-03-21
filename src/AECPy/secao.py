'''
Módulo para definição da classe Secao para o AECPy
'''

import math
from material import Material

class Secao():

    '''
    Seção transversal para um elemento de barra prismático
    usado na análise estrutural por AECPy

    Na denominação dos atributos, são considerados:
    eixo local 2 - eixo vertical
    eixo local 3 - eixo horizonal
    Os eixos 2 e 3 são eixos principais centrais de inércia da seção
    
    Atributes
    ----------
    mat: Material
        Material que forma o elemento estrutural
    A: float
        Área da seção transveral
    I2: float
        Momento de inércia em relação ao eixo local 2
    I3: float
        Momento de inércia em relação ao eixo local 3
    J: float
        Constante de torção pura (constante de St. Venant)
    AE: float
        Rigidez axial da seção
    EI2: float
        Rigidez da seção à flexão em torno do eixo local 2
    EI3: float
        Rigidez da seção à flexão em torno do eixo local 3 
    GJ: float
        Rigidez à torção pura da seção
    nome: str
        Nome do material
    peso_unitario
    r2
    r3
    
    '''

    def __init__(self,mat,A,I3=0,I2=0,J=0,W3=0,W2=0,nome=''):
        '''
        mat: Material
            Material que forma o elemento estrutural
        A: float
            Área da seção transveral
        I2: float (optional)
            Momento de inércia em relação ao eixo local 2. (deault é 0, aplicável à análise de treliças)
        I3: float (optional)
            Momento de inércia em relação ao eixo local 3 (default é 0, não aplicável para análise de pórticos espaciais)
        J: float (optional)
            Constante de torção pura (constante de St. Venant) (default é 0, não aplicável à análise de pórticos espaciais)
        W3: float(optional)
            Módulo de resistência (elástico) para flexão em relação ao eixo 3 (default 0, não aplicável para análise com variaçao térmica na direção do eixo 2)
        W2: float (optional)
            Módulo de resistência (elástico) para flexão em relação ao eixo 2 (default 0, não aplicável para análise com variaçao térmica na direção do eixo 3)
        nome: str (optional)
            Nome da seção transveral (default "")

        Raises
        ------
        ValueError
            Se qualquer propriedade geométrica for negativa
        '''
        
        if not isinstance(mat,Material):
            raise TypeError("mat deve ser do tipo Material")
        if any([ prop<0 for prop in [A,I2,I3,J]]):
                raise ValueError('As propriedades geométricas não devem ser negativas')

        self.mat = mat   #material 
        self.A = A      #Área da seção
        self.I2 = I2    #Momento de inérica da seção em relação ao eixo local 2
        self.I3 = I3    #Momento de inérica da seção em relação ao eixo local 3
        self.J = J      #Constante de torção pura (St. Venant)            
        self.W2 = W2    #Módulo de resistência elástico em relação ao eixo 2
        self.W3 = W3    #Módulo de resistência elástico em relação ao eixo 3
        self.nome = nome

        self.EA = mat.E*self.A     #Rigidez da seção ao alongamento/encurtamento axial
        self.EI2 = mat.E*self.I2   #Ridigez da seção à flexão em torno do eixo local 2 
        self.EI3 = mat.E*self.I3   #Ridigez da seção à flexão em torno do eixo local 3
        self.GJ = mat.G*self.J     #Rigidez da seção à torção pura (St. Venant)


    @property
    def peso_unitario(self):
        '''Peso por unidade de comprimento (da barra)'''
        return self.mat.pe * self.A
    
    @property
    def r2(self):
        '''Raio de giração em relação ao eixo 2'''
        return math.sqrt(self.I2/self.A)
    
    @property
    def r3(self):
        '''Raio de giração em relação ao eixo 3'''
        return math.sqrt(self.I3/self.A)


    def __repr__(self):
        return f'Secao(mat={self.mat.nome}, A ={self.A}, I2={self.I2}, I3={self.I3},J ={self.J} )'
    
    def __str__(self):
        txt  = f'secao: {self.nome} \n'
        txt += f'A ={self.A} \n'
        txt += f'I2={self.I2} \n'
        txt += f'I3={self.I3} \n'
        txt += f'J ={self.J} \n'
        return txt
    

class SecaoRetangular(Secao):
    '''Seção retangular para análise estrutural pelo AECPy'''

    def __init__(self,mat,b,h,nome=''):
        '''
        Parameters
        ----------
        mat: Material
            Material da seção transveral
        b: float
            Largura da seção. Comprimento lado paralelo ao eixo 2.
        h: float
            Altura da seção. Comprimento do paralelo ao eixo 3.
        nome: str (optional)
            Nome da seção transveral (default "")
        '''

        A = b*h
        I2 = b**3*h/12
        I3 = b*h**3/12
        #constante de torção pura
        a = max(h,b)/2
        c = min(h,b)/2
        J  = (a*c**3)*(16/3 - 3.36*(c/a)*(1-c**4/(12*a**4)) )
        
        W2 = I2 / (b/2)
        W3 = I3 / (h/2)
        self.h = h
        self.b = b
        super().__init__(mat,A=A,I3=I3,I2=I2,J=J,W2=W2,W3=W3,nome=nome)


