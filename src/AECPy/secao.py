"""
Módulo para definição das classes de material e seção transversal para o AECPy
"""


class Material:
    """Material elástico linear isotrópico para análise estrutural.

    Parameters
    ----------
    E : float
        Módulo de Young.
    G : float, optional
        Módulo de cisalhamento (default 0).
    pe : float, optional
        Peso específico (default 0).
    cdt : float, optional
        Coeficiente de dilatação térmica (default 0).
    nome : str, optional
        Nome do material (default "").

    Examples
    --------
    >>> Material(E=200e9, G=77e9, pe=78500, cdt=1.2e-5, nome='Aço')
    """

    def __init__(self, E, G=0.0, pe=0.0, cdt=0.0, nome=""):
        if any(prop < 0 for prop in [E, G, pe, cdt]):
            raise ValueError("As propriedades do material não devem ser negativas")

        self.E = E          # módulo de Young
        self.G = G          # módulo de cisalhamento
        self.pe = pe        # peso específico
        self.cdt = cdt      # coeficiente de dilatação térmica
        self.nome = nome

    def __repr__(self):
        return f"Material({self.E}, {self.G}, {self.pe}, {self.cdt}, {self.nome!r})"

    def __str__(self):
        txt = "material" if self.nome == "" else self.nome
        txt += f": E={self.E}, G={self.G}"
        txt += f", peso específico={self.pe}" 
        txt += f", coef. dilatação térmica={self.cdt}"
        return txt


class SecaoBase:
    """Seção transversal definida diretamente pelas rigidezes e coeficientes térmicos.

    Útil quando as rigidezes são conhecidas diretamente (p. ex., perfis de catálogo),
    sem necessidade de fornecer material ou propriedades geométricas.

    Parameters
    ----------
    EA : float
        Rigidez axial.
    EI3 : float, optional
        Rigidez à flexão em torno do eixo local 3 (default 0).
    EI2 : float, optional
        Rigidez à flexão em torno do eixo local 2 (default 0).
    GJ : float, optional
        Rigidez à torção pura, St. Venant (default 0).
    peso_unitario : float, optional
        Peso próprio por unidade de comprimento (default 0).
    cdtEA : float, optional
        Produto ``alpha * EA`` para efeitos térmicos axiais (default 0).
    cdtEI2 : float, optional
        Produto ``alpha * EI2`` para efeitos térmicos na flexão em torno do eixo 2 (default 0).
    cdtEI3 : float, optional
        Produto ``alpha * EI3`` para efeitos térmicos na flexão em torno do eixo 3 (default 0).
    nome : str, optional
        Nome da seção (default "").

    Examples
    --------
    >>> SecaoBase(EA=4e9, EI3=1.3e7, EI2=3.3e6, GJ=4.3e5)
    """

    def __init__(self, EA, EI3=0.0, EI2=0.0, GJ=0.0, peso_unitario=0.0,
                 cdtEA=0.0, cdtEI2=0.0, cdtEI3=0.0, nome=""):
        if any(v < 0 for v in [EA, EI2, EI3, GJ, peso_unitario, cdtEA, cdtEI2, cdtEI3]):
            raise ValueError("Os atributos de SecaoBase não devem ser negativos")
        self.EA = EA
        self.EI2 = EI2
        self.EI3 = EI3
        self.GJ = GJ
        self.peso_unitario = peso_unitario
        self.cdtEA = cdtEA
        self.cdtEI2 = cdtEI2
        self.cdtEI3 = cdtEI3
        self.nome = nome

    def __repr__(self):
        return (f"SecaoBase(EA={self.EA}, EI3={self.EI3}, EI2={self.EI2}, "
                f"GJ={self.GJ})")

    def __str__(self):
        txt = f"secao: {self.nome}\n"
        txt += f"EA={self.EA:.3e}  EI3={self.EI3:.3e}  EI2={self.EI2:.3e}  GJ={self.GJ:.3e}\n"
        txt += f"cdtEA={self.cdtEA:.3e}  cdtEI3={self.cdtEI3:.3e}  cdtEI2={self.cdtEI2:.3e}  peso={self.peso_unitario:.3e}\n"
        return txt


class Secao(SecaoBase):
    """Seção transversal de barra prismática definida por material e propriedades geométricas.

    Calcula as rigidezes e coeficientes térmicos a partir de ``mat``, ``A``,
    ``I2``, ``I3`` e ``J``. Eixos locais: eixo 2 vertical, eixo 3 horizontal
    (eixos principais centrais de inércia).

    Parameters
    ----------
    mat : Material
        Material do elemento.
    A : float
        Área da seção transversal.
    I3 : float, optional
        Momento de inércia em relação ao eixo local 3 (default 0).
    I2 : float, optional
        Momento de inércia em relação ao eixo local 2 (default 0).
    J : float, optional
        Constante de torção pura, St. Venant (default 0).
    nome : str, optional
        Nome da seção (default "").

    Examples
    --------
    >>> Secao(mat, A=0.02, I3=6.7e-5, I2=1.7e-5, J=5.6e-6)
    """

    def __init__(self, mat, A, I3=0.0, I2=0.0, J=0.0, nome=""):
        if not isinstance(mat, Material):
            raise TypeError("mat deve ser do tipo Material")

        self.mat = mat
        self.A = A
        self.I2 = I2
        self.I3 = I3
        self.J = J

        EA  = mat.E * A
        EI2 = mat.E * I2
        EI3 = mat.E * I3
        GJ  = mat.G * J

        super().__init__(
            EA=EA, EI2=EI2, EI3=EI3, GJ=GJ,
            peso_unitario=mat.pe * A,
            cdtEA=mat.cdt * EA,
            cdtEI2=mat.cdt * EI2,
            cdtEI3=mat.cdt * EI3,
            nome=nome,
        )

    def __repr__(self):
        return (f"Secao(mat={self.mat.nome}, A={self.A}, "
                f"I2={self.I2}, I3={self.I3}, J={self.J})")

    def __str__(self):
        txt = super().__str__()
        txt += str(self.mat) + "\n"
        txt += f"A={self.A:.4g}  I2={self.I2:.3e}  I3={self.I3:.3e}  J={self.J:.3e}\n"
        return txt


class SecaoRetangular(Secao):
    """Seção retangular de barra prismática.

    Calcula ``A``, ``I2``, ``I3`` e ``J`` (St. Venant) a partir de ``b`` e ``h``.

    Parameters
    ----------
    mat : Material
        Material do elemento.
    b : float
        Largura da seção (paralelo ao eixo local 2).
    h : float
        Altura da seção (paralelo ao eixo local 3).
    nome : str, optional
        Nome da seção (default "").

    Examples
    --------
    >>> SecaoRetangular(mat, b=0.1, h=0.2)
    """

    def __init__(self, mat, b, h, nome=""):
        A  = b * h
        I2 = b**3 * h / 12
        I3 = b * h**3 / 12
        # constante de torção pura (St. Venant)
        a = max(h, b) / 2
        c = min(h, b) / 2
        J = (a * c**3) * (16/3 - 3.36 * (c/a) * (1 - c**4 / (12 * a**4)))

        self.b = b
        self.h = h
        super().__init__(mat, A=A, I3=I3, I2=I2, J=J, nome=nome)

    def __str__(self):
        txt = super().__str__()
        txt += f"b={self.b:.4g}  h={self.h:.4g}\n"
        return txt
