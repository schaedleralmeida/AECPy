"""
Módulo para definição da classe material para o AECPy
"""


class Material:
    """
    Material elástico linear para análise estrutural por AECPy

    Atributes
    ----------
    E: float
        Módulo de Young
    G: float
        Módulo de cisalhamento
    cp: float
        Coeficiente de Poisson
    pe: float
        Peso específico
    cdt: float
        Coeficiente de dilatação térmica
    nome: str
        Nome do material
    """

    def __init__(self, E, cp, pe=0.0, cdt=0.0, nome=""):
        """
        Inicialização de uma instância da classe Material

        Parametres
        ----------
        E: float
            Módulo de Young
        cp: float
            Coeficiente de Poisson
        pe: float (optional)
            Peso específico (default 0)
        cdt: float (optional)
            Coeficiente de dilatação térmica (default 0)
        nome: str (optional)
            Nome do material (default "")

        Raises
        ------
        ValueError
            Se qualquer propriedade do material for negativa
        """

        if any([prop < 0 for prop in [E, cp, pe, cdt]]):
            raise ValueError(
                "As propriedades do material não devem ser negativas"
            )

        self.E = E  # módulo de Young (de elasticidade)
        self.cp = cp  # coeficiente de poisson
        self.G = E / (
            2 * (1 + cp)
        )  # módulo de cisalhamento (de elasticidade transversal)
        self.pe = pe  # peso específico
        self.cdt = cdt  # Coefciente de dilatação térmica
        self.nome = nome  # nome do material

    def __repr__(self):
        return f"material({self.E},{self.cp},{self.pe},{self.cdt},{self.nome})"

    def __str__(self):
        txt = "material" if self.nome == "" else self.nome
        txt += f": E={self.E}, Poisson={self.cp}, G={self.G}"
        txt += f", peso específico={self.pe}" if self.pe != 0 else ""
        txt += f", coef. dilatação térmica={self.cdt}" if self.cdt != 0 else ""
        return txt
