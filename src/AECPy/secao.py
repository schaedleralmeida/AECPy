"""
Módulo para definição da classe Secao para o AECPy
"""


class Material:
    """Material elástico linear para análise estrutural por AECPy.

    Atributos
    ---------
    E : float
        Módulo de Young.
    G : float
        Módulo de cisalhamento.
    pe : float
        Peso específico.
    cdt : float
        Coeficiente de dilatação térmica.
    nome : str
        Nome do material.
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
    """Classe base de seção transversal para elementos estruturais.

    Armazena diretamente as rigidezes e os coeficientes térmicos da seção,
    sem dependência do material ou das propriedades geométricas brutas.
    Útil quando as rigidezes são conhecidas diretamente (p. ex., perfis de catálogo).

    Atributos
    ---------
    EA : float
        Rigidez axial.
    EI2 : float
        Rigidez à flexão em torno do eixo local 2.
    EI3 : float
        Rigidez à flexão em torno do eixo local 3.
    GJ : float
        Rigidez à torção pura (St. Venant).
    peso_unitario : float
        Peso próprio por unidade de comprimento.
    cdtEA : float
        Coeficiente de dilatação térmica axial (``= alpha * EA``).
    cdtEI2 : float
        Coeficiente de dilatação térmica à flexão no eixo 2 (``= alpha * EI2``).
    cdtEI3 : float
        Coeficiente de dilatação térmica à flexão no eixo 3 (``= alpha * EI3``).
    nome : str
        Nome da seção.
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
    """Seção transversal para um elemento de barra prismático.

    Calcula as rigidezes e os coeficientes térmicos a partir das propriedades
    do material e das dimensões geométricas da seção.

    Na denominação dos eixos locais:
    - eixo 2: eixo vertical (principal central de inércia).
    - eixo 3: eixo horizontal (principal central de inércia).

    Atributos adicionais (além dos de SecaoBase)
    --------------------------------------------
    mat : Material
        Material do elemento.
    A : float
        Área da seção transversal.
    I2 : float
        Momento de inércia em relação ao eixo local 2.
    I3 : float
        Momento de inércia em relação ao eixo local 3.
    J : float
        Constante de torção pura (St. Venant).
    """

    def __init__(self, mat, A, I3=0.0, I2=0.0, J=0.0, nome=""):
        """
        Parameters
        ----------
        mat : Material
            Material que forma o elemento estrutural.
        A : float
            Área da seção transversal.
        I3 : float, optional
            Momento de inércia em relação ao eixo local 3 (default 0).
        I2 : float, optional
            Momento de inércia em relação ao eixo local 2 (default 0).
        J : float, optional
            Constante de torção pura (St. Venant) (default 0).
        nome : str, optional
            Nome da seção transversal (default "").

        Raises
        ------
        TypeError
            Se ``mat`` não for do tipo Material.
        ValueError
            Se qualquer propriedade geométrica for negativa.
        """
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
    """Seção retangular para análise estrutural pelo AECPy."""

    def __init__(self, mat, b, h, nome=""):
        """
        Parameters
        ----------
        mat : Material
            Material da seção transversal.
        b : float
            Largura da seção (paralelo ao eixo 2).
        h : float
            Altura da seção (paralelo ao eixo 3).
        nome : str, optional
            Nome da seção transversal (default "").
        """
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
