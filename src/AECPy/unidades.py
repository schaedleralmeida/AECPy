"""
Módulo para conversão de unidades de medidas.

O objetivo principal é a transformação de quantiades entre
unidades defindidas pelo usário e unidades definidas como padrão
para os cálculos em um programa.

Classes:

    Quantidade: se destina à conversão de unidades de medidas
    para uma determinada quantidade (ex.: comrpimento, força, momento, etc.)
    Em cada quantidade estão definidas unidades de medidas,
    os fatores de transformação entre estas unidades e
    a unidade padrão usada nos cálculos

    Conversor: se destina à conversão de unidades para um conjunto
    de quantidades usadas nos cálculos. As unidades padrão de cada quantidade,
    são compatíveis com as unidades padrão de comprimento e força
    para a realização dos cálculos

funções:

    unidades_derivadas_expoente(unidades,expoente)

    unidades_derivadas_prefixos(unidade,fator,prefixos)

    unidades_derivadas_compor(unidades1,unidades2,operacao='multiplicar')

"""

__all__ = ["Quantidade", "Conversor"]
__version__ = "0.1"
__author__ = "Felipe Schaedler de Almeida"


class Quantidade:
    """
    Classe para definir a conversão entre unidades de uma quantidade

    Atributos
    ---------

    Propriedades
    ------------
    nome: str
        Nome da quantidade definida (ex. `força`)
    padrao: str
        Unidade padrão usada na conversão (ex. `kN` para força)
    unidades: dict
        Unidades contempladas para a quantidade
        key: é o nome da unidade
        value: é o fator de conversao para a unidade básica do SI
        exemplo: {'N':1, 'kN':1000, 'kgf':10} para força

    Exemplo
    -------
    quandidade de comprimento, tendo milímetros ('mm') como unidade padrão
    e considerando as unidades metros ('m') e centímetros ('cm') para convresão:
    >>> q = Quantidade('comprimento','mm',{'m':1, 'cm':0.01, 'mm':0.001})

    conversão de 200 cm para a unidade padrão ('mm'):
    >>> q = Quantidade('comprimento','mm',{'m':1, 'cm':0.01, 'mm':0.001})
    >>> 200 q.em('cm')
    2000.0

    conversão da quandidade na unidade padrão ('mm') para a unidade definida ('m'):
    >>> q = Quantidade('comprimento','mm',{'m':1, 'cm':0.01, 'mm':0.001})
    >>> 2000.0 * q.para('m')
    2.0

    criação de uma string com o valor e a unidade padrão:
    >>> q = Quantidade('comprimento','mm',{'m':1, 'cm':0.01, 'mm':0.001})
    >>> q.para_str(2000.0)
    '2000.0 mm'

    criação de uma string com valor e a unidade em 'cm'
    >>> q = Quantidade('comprimento','mm',{'m':1, 'cm':0.01, 'mm':0.001})
    >>> p.para_str(2000.0,'cm')
    '200 cm'

    obter a unidade padrão da quandidade
    >>> q = Quantidade('comprimento','mm',{'m':1, 'cm':0.01, 'mm':0.001})
    >>> q.padrao
    'mm'

    obter as unidades definidas para a conversão
    >>> q = Quantidade('comprimento','mm',{'m':1, 'cm':0.01, 'mm':0.001})
    >>> q.unidades
    {'m':1, 'cm':0.01, 'mm':0.001}
    """

    def __init__(self, nome, unidade_padrao, unidades):
        """
        Inicia uma instânica da classe quantidade

        Parâmetros
        ----------
        nome: str
            Nome da quantidade
        unidade_padrao: str
            Unidade padrão para coversão (ex. `kN`)
        unidades: dict
            dicionário com unidades, onde
            key: é o nome da unidade
            value: é o fator de conversao para a unidade básica do SI
        """

        self.__nome = nome
        self.__unidades = unidades
        self.definir_padrao(unidade_padrao)

    @property
    def nome(self):
        """Nome da quantidade"""
        return self.__nome

    @property
    def unidades(self):
        """Unidades definidas para a quantidade"""
        return self.__unidades

    @property
    def padrao(self):
        """Unidade padrão para a quantidade"""
        return self.__padrao

    def __verificar(self, uni):
        """Verifica se a unidade `uni` está definida para a quantidade"""
        if uni not in self.unidades:
            raise ValueError(
                f"A unidade {uni} não está na lista de unidades {list(self.unidades.keys())}"
            )

    def definir_padrao(self, unidade_padrao):
        """Define a `unidade_padrao` como a unidade padrão da quantidade"""
        self.__verificar(unidade_padrao)
        self.__padrao = unidade_padrao
        self.__conv = dict()
        for uni in self.unidades:
            self.__conv[uni] = (
                self.unidades[uni] / self.unidades[unidade_padrao]
            )

    def para_padrao(self, val, uni):
        """Converte o valor `val` da unidade `uni` para a unidade padrão"""
        self.__verificar(uni)
        return val * self.__conv[uni]

    def de_padrao(self, val, uni):
        """Converte o valor `val` da unidade padrão para a unidade `uni`"""
        self.__verificar(uni)
        return val / self.__conv[uni]

    def de_para(self, val, de_uni, para_uni):
        """Converte o valor `val` da unidade `de_uni` para a unidade `para_uni`"""
        val_uni_padrao = self.para_padrao(val, de_uni)
        return self.de_padrao(val_uni_padrao, para_uni)

    def para_str(self, val, uni=None):
        """
        Converte o valor `val` da unidade padrão para a unidade `uni`
        e cria uma string com o valor e a unidade.
        Se `uni` é None, adota a unidade padrão da quantidade
        """
        if uni is None:
            valor = val
            unidade = self.padrao
        else:
            valor = self.de_padrao(val, uni)
            unidade = uni
        return str(valor) + " " + unidade

    def em(self, uni):
        """
        Retorna o fator de conversao da unidade `uni` para a unidade padrão

        Parâmetro
        ---------
        uni: str
            Unidade em que a quantidade é definida

        Exemplo
        -------
        >>> q = Quantidade('comprimento','m',{'m':1, 'cm':0.01, 'mm':0.001})
        >>> 200 * q.em('mm') : converte 200 mm para a unidade padrão de comprimento (2.0)
        """
        self.__verificar(uni)
        return self.__conv[uni]

    def para(self, uni):
        """
        Retorna o fator de conversao da unidade padrão para a unidade `uni`

        Parâmetro
        ---------
        uni: str
            Unidade em que a quantidade é definida

        Exemplo
        -------
        >>> q = Quantidade('comprimento','m',{'m':1, 'cm':0.01, 'mm':0.001})
        >>> 2 * q.para('mm') : converte 2 em unidade padrão de comprimento ('m') para mm (2000)
        """
        self.__verificar(uni)
        return 1 / self.__conv[uni]

    def __repr__(self) -> str:
        return f"quantidade({self.nome}, {self.padrao}, {self.unidades})"

    def __str__(self) -> str:
        txt = f"nome: {self.nome} \n"
        txt += f"padrao: {self.padrao} \n"
        txt += f"unidades: {self.unidades.__repr__()}"
        return txt


def unidades_derivadas_expoente(unidades, expoente):
    """
    Cria um dicionário de unidades derivadas pelo expoente

    Parâmetros
    ----------
    unidades: dict
        dicionário com unidades, onde
        key: é o nome da variável
        value: é o fator de conversao para a unidade básica do SI
    expoente: int
        Define o expoente da unidade derivada

    Returns
    -------
    uni: dict
        dicionário com unidades, onde
        key: é o nome da variável
        value: é o fator de conversao para a unidade básica do SI
    """

    exp = str(expoente)
    uni = dict()
    for u, f in unidades.items():
        uni[u + exp] = f**expoente
    return uni


def unidades_derivadas_compor(unidades1, unidades2, operacao="multiplicar"):
    """
    Cria um dicionário de unidades derivadas pela composição de duas variáveis

    Parâmetros
    ----------
    unidades1, unidades2: dict
        dicionários com unidades, onde
        key: é o nome da variável
        value: é o fator de conversao para a unidade básica do SI
    operacao: str (opcional, default é multiplicar)
        Define como as unidades vão compor a unidade derivada.
        São aceitos os valores:
            `multiplicar` -> unidades1 * unidades2
            `dividir` -> unidades1/ unidades2

    Returns
    -------
    uni: dict
        dicionário com unidades, onde
        key: é o nome da variável
        value: é o fator de conversao para a unidade básica do SI
    """

    uni = dict()
    for u1, f1 in unidades1.items():
        for u2, f2 in unidades2.items():
            if operacao == "multiplicar":
                uni[u1 + u2] = f1 * f2
            elif operacao == "dividir":
                uni[u1 + "/" + u2] = f1 / f2
    return uni


def unidades_derivadas_prefixos(unidade, fator, prefixos):
    """
    Cria um dicionário de unidades derivadas em termos do prefixo

    Parâmetros
    ----------
    unidade: str
        Unidade básica usada para crirar as unidades derivadas
    fator: float
        Fator multiplicador da unidade básica em relação à unidade padrao do SI
    prefixos: list
        Lista de prefixos usados para criar as unidades derivdas
        São válidos os prefixos k, M e G


    Returns
    -------
    uni: dict
        dicionário com unidades, onde
        key: é o nome da variável
        value: é o fator de conversao para a unidade básica do SI
    """
    dp = {"k": 10**3, "M": 10**6, "G": 10**9}
    uni = {unidade: fator}
    for pref in prefixos:
        if pref in dp:
            uni[pref + unidade] = dp[pref] * fator
    return uni


class Conversor:
    """
    Classe para transformação de unidadades de várias quantidades
    usadas em um mesmo problema.

    Propriedades
    ------------
    quantidades: list
        Lista das quantidade contemploadas no conversor
    unidades_padrao: dict
        Dicionario com a unidade padrão para cada quantidade considerada no conversor
    comprimento, forca, momento, area , volume, tensao, peso_especifico: quantidade
        acesso direto às quantidades disponíveis para conversão de unidades

    Exemplo
    -------
    conversor, tendo mm como undiade padrão de comprimento e 'N' como unidade padrão de força:
    >>> c = Conversor(up_comprimento='mm', up_forca='N')

    conversão de 200 cm para a unidade padrão ('mm'):
    >>> c = Conversor(up_comprimento='mm', up_forca='N')
    >>> 200 * c.em('cm')
    2000.0
    >>> 200 * c.comprimento.em('cm')
    2000.0

    conversão de 15  kNm para a unidade padrão ('Nmm'):
    >>> c = Conversor(up_comprimento='mm', up_forca='N')
    >>> 15 * c.em(kNm)
    1.5e7
    >>> 15 * c.momento.em('kNm)
    1.5e7

    conversão de momento na unidade padrão ('Nmm') para 'kNcm':
    >>> c = Conversor(up_comprimento='mm', up_forca='N')
    >>> 15.e6 * c.para('kNcm')
    1500.0
    >>> 15.e6 *c.momento.para('kNcm')
    1500.0

    criação de uma string com o valor e a unidade padrão:
    >>> c = Conversor(up_comprimento='mm', up_forca='N')
    >>> c.para_str(2000.0,'comprimento')
    '2000.0 mm'
    >>> c.comprimento.para_str(2000.0)
    '2000.0 mm'

    criação de uma string com valor e a unidade em 'cm'
    >>> c = Conversor(up_comprimento='mm', up_forca='N')
    >>> c.para_str(2000.0,'cm')
    '200.0 cm'
    >>> c.comprimento.para_str(2000.0,'cm')
    '200.0 cm'

    lista das unidades padraão das quantidades disponíveis para convresão
    >>> c = Conversor(up_comprimento='mm', up_forca='N')
    >>> c.unidades_padrao
    {'momento': 'Nmm',
    'força': 'N',
    'comprimento': 'mm',
    'área': 'mm2',
    'volume': 'mm3',
    'tensão': 'N/mm2',
    'peso específico': 'N/mm3'}
    """

    def __init__(self, up_comprimento="m", up_forca="N"):
        self.__quantidades = list()

        # comprimento
        self.__comprimento = Quantidade(
            "comprimento", up_comprimento, {"m": 1, "mm": 1.0e-3, "cm": 1.0e-2}
        )
        self.__quantidades.append(self.__comprimento)

        # força
        self.__forca = Quantidade(
            "força", up_forca, {"N": 1, "kN": 1.0e3, "kgf": 1.0e1}
        )
        self.__quantidades.append(self.__forca)

        # momento
        unidades = unidades_derivadas_compor(
            self.__forca.unidades, self.__comprimento.unidades
        )
        up = up_forca + up_comprimento
        self.__momento = Quantidade("momento", up, unidades)
        self.__quantidades.append(self.__momento)

        # área
        unidades = unidades_derivadas_expoente(self.__comprimento.unidades, 2)
        up = up_comprimento + "2"
        self.__area = Quantidade("área", up, unidades)
        self.__quantidades.append(self.__area)

        # volume
        unidades = unidades_derivadas_expoente(self.__comprimento.unidades, 3)
        up = up_comprimento + "3"
        self.__volume = Quantidade("volume", up, unidades)
        self.__quantidades.append(self.__volume)

        # tensão
        unidades = unidades_derivadas_compor(
            self.__forca.unidades, self.__area.unidades, "dividir"
        )
        unidades_adicionais = unidades_derivadas_prefixos(
            "Pa", 1, ["k", "M", "G"]
        )
        up = up_forca + "/" + up_comprimento + "2"
        self.__tensao = Quantidade(
            "tensão", up, unidades | unidades_adicionais
        )
        self.__quantidades.append(self.__tensao)

        # peso específico
        unidades = unidades_derivadas_compor(
            self.__forca.unidades, self.__volume.unidades, "dividir"
        )
        up = up_forca + "/" + up_comprimento + "3"
        self.__peso_especifico = Quantidade("peso específico", up, unidades)
        self.__quantidades.append(self.__peso_especifico)

        return

    @property
    def quantidades(self):
        """Quantidades definidas para conversão"""
        return self.__quantidades

    @property
    def unidades_padrao(self):
        """Unidades padrão para conversão das quantidades definidas"""
        up = dict()
        for q in self.quantidades:
            up[q.nome] = q.padrao
        return up

    @property
    def comprimento(self):
        """quantidade comprimento"""
        return self.__comprimento

    @property
    def forca(self):
        """quantidade força"""
        return self.__forca

    @property
    def area(self):
        """quantidade área"""
        return self.__area

    @property
    def volume(self):
        """quantidade volume"""
        return self.__volume

    @property
    def momento(self):
        """quantidade momento"""
        return self.__momento

    @property
    def tensao(self):
        """quantidade tensão"""
        return self.__tensao

    @property
    def peso_especifico(self):
        """quantidade peso específico"""
        return self.__peso_especifico

    def get_quantidade(self, uni):
        """retorna a quantidade que contém a unidade `uni` e levanta ValueErro caso não encontre"""
        for i, q in enumerate(self.quantidades):
            if uni in q.unidades:
                if i > 0:
                    self.__quantidades.insert(0, self.__quantidades.pop(i))
                return q
        raise ValueError(f"A unidade {uni} não está na lista de unidades")

    def em(self, uni):
        """
        Retorna o fator de conversao da unidade `uni` para a unidade padrão

        Parâmetro
        ---------
        uni: str
            Unidade em que a quantidade é definida

        Exemplo
        -------
        >>> conv = Conversor()
        >>> 200 * conv.em('mm') : converte 200 mm para a unidade padrão de comprimento
        """
        q = self.get_quantidade(uni)
        return q.em(uni)

    def para(self, uni):
        """
        Retorna o fator de conversao da unidade padrão para a unidade `uni`

        Parâmetro
        ---------
        uni: str
            Unidade em que a quantidade é definida

        Exemplo
        -------
        >>> conv = Conversor()
        >>> 200 * conv.para('mm') : converte 200 em unidade padrão de comprimento para mm
        """
        q = self.get_quantidade(uni)
        return q.para(uni)

    def para_str(self, val, uni_nome):
        """
        Converte o valor em uma string com a unidade

        Parâmetros
        ----------
        val: float
            Valor da quantidade na unidade padrão
        uni_nome: str
            Define a unidade ou tipo de quantidade.
            Quando a unidade é informada, o valor `val` é convertido
            Quando a quantidade é informada, a unidade padrão é usada para formar a string
        """
        for q in self.quantidades:
            if q.nome == uni_nome:
                return q.para_str(val)
        q = self.get_quantidade(uni_nome)
        return q.para_str(val, uni_nome)
