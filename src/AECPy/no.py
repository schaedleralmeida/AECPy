"""
Módulo para definição dos nós na análise estrutural
pelo Método da Rigidez Direta (MRD) no AECPy
"""

from numbers import Integral, Real

import numpy as np


class NoBase:
    """Classe base de nó para modelos estruturais via Método da Rigidez Direta.

    Cada subclasse define o tipo de modelo (PP, PE, TP, TE, GR), incluindo
    os graus de liberdade e componentes de força disponíveis no nó.
    """

    @staticmethod
    def _eh_numero(valor) -> bool:
        return isinstance(valor, Real) and not isinstance(valor, bool)

    @staticmethod
    def _eh_inteiro_natural(valor) -> bool:
        return (
            isinstance(valor, Integral)
            and not isinstance(valor, bool)
            and valor >= 0
        )

    # atributos do tipo estrutural (definidos em subclasses)
    tipo = ""
    ndim = 0
    ngdl = 0
    eixos_locais = ()
    eixos_globais = ()
    gdls_globais = ()
    gdls_locais = ()
    forcas_globais = ()


    def __init__(self, coor) -> None:
        if self.ndim <= 0 or self.ngdl <= 0:
            raise TypeError(
                "A classe NoBase é base. Use uma subclasse como NoPE, NoPP, NoTE, NoTP ou NoGR."
            )
        if len(coor) != self.ndim:
            raise ValueError("Número incorreto de coordenadas")
        if not all([self._eh_numero(c) for c in coor]):
            raise TypeError("As coordenadas devem ser números")

        # atributos de instância
        self.coor = np.array(coor)  # coordenadas da posição do nó
        self._igdl = ()  # índice global dos gdl do nó
        self.__deslocamentos_nulos = ()  # nomes dos gdl com apoios (deslocamentos nulos)
        self.__deslocamentos_prescritos = {}
        self.__apoio_elastico = {}
        self.__carga = {}

    @property
    def x(self):
        """Coordenada global x do nó."""
        return self.coor[0]

    @property
    def y(self):
        """Coordenada global y (quando o modelo possui eixo y)."""
        return self.coor[self.eixos_globais.index("y")]

    @property
    def z(self):
        """Coordenada global z (quando o modelo possui eixo z)."""
        return self.coor[self.eixos_globais.index("z")]

    @property
    def igdl(self):
        """Índices globais dos graus de liberdade do nó no sistema montado."""
        return list(self._igdl)

    @igdl.setter
    def igdl(self, ig):
        if len(ig) != self.ngdl:
            raise ValueError("O número de gdl está incorreto")
        if not all(
            [
                self._eh_inteiro_natural(i)
                for i in ig
            ]
        ):
            raise ValueError(
                f"Os índices globais dos gdl do nó devem ser números naturais: {ig}"
            )
        self._igdl = tuple(ig)

    def idl(self, deslocamento):
        """Retorna o índice local de um deslocamento (ex.: ux, uz, ry)."""
        if deslocamento not in self.gdls_globais:
            raise ValueError(f"Deslocamento inválido: {deslocamento}")
        return self.gdls_globais.index(deslocamento)

    def idg(self, deslocamento):
        """Retorna o índice global do deslocamento para o nó atual."""
        return self._igdl[self.idl(deslocamento)]

    def ifl(self, forca):
        """Retorna o índice local de uma componente de força/momento."""
        if forca not in self.forcas_globais:
            raise ValueError(f"Força inválida: {forca}")
        return self.forcas_globais.index(forca)

    def ifg(self, forca):
        """Retorna o índice global da componente de força/momento no nó."""
        return self._igdl[self.ifl(forca)]

    @property
    def deslocamentos_nulos(self):
        """Deslocamentos com valor nulo por apoio (ex.: ['ux', 'uz'])."""
        return list(self.__deslocamentos_nulos)

    @deslocamentos_nulos.setter
    def deslocamentos_nulos(self, valores):
        """Define deslocamentos nulos por gdl.

        Exemplo: ['ux', 'uz']; use None para limpar os apoios.
        """
        if valores is None:
            self.__deslocamentos_nulos = ()
        else:
            # Validar que todos os valores estão em gdls_globais
            if not isinstance(valores, (list, tuple, set)):
                raise TypeError("deslocamentos_nulos deve ser uma lista, tupla ou conjunto")
            for valor in valores:
                self.idl(valor)
            self.__deslocamentos_nulos = tuple(valores)
        
        # Verificar conflitos de condições de contorno
        self._check_cc()

    @property
    def deslocamentos_prescritos(self):
        """Dicionário de deslocamentos prescritos por gdl.

        Exemplo: {'ux': 0.002, 'ry': -1e-3}
        """
        return dict(self.__deslocamentos_prescritos)

    @deslocamentos_prescritos.setter
    def deslocamentos_prescritos(self, valores):
        """Define deslocamentos prescritos por gdl.

        Exemplo: {'ux': 0.002, 'ry': -1e-3}; use None para limpar.
        """
        if valores is None:
            self.__deslocamentos_prescritos = {}
            self._check_cc()
            return

        if not isinstance(valores, dict):
            raise TypeError("deslocamentos_prescritos deve ser um dicionário")

        for gdl, valor in valores.items():
            self.idl(gdl)
            if not self._eh_numero(valor):
                raise TypeError(
                    f"O valor do deslocamento prescrito em '{gdl}' deve ser numérico"
                )

        self.__deslocamentos_prescritos = dict(valores)

        # Verificar conflitos de condições de contorno
        self._check_cc()

    @property
    def apoio_elastico(self):
        """Dicionário de apoios elásticos por gdl.

        Cada entrada associa um gdl à sua rigidez de mola equivalente.
        """
        return dict(self.__apoio_elastico)

    @apoio_elastico.setter
    def apoio_elastico(self, valores):
        """Define rigidezes de apoio elástico por gdl.

        Exemplo: {'ux': 1.0e5, 'ry': 2.5e6}; use None para limpar.
        """
        if valores is None:
            self.__apoio_elastico = {}
            self._check_cc()
            return

        if not isinstance(valores, dict):
            raise TypeError("apoio_elastico deve ser um dicionário")

        for gdl, valor in valores.items():
            self.idl(gdl)
            if not self._eh_numero(valor):
                raise TypeError(
                    f"O valor do apoio elástico em '{gdl}' deve ser numérico"
                )

        self.__apoio_elastico = dict(valores)

        # Verificar conflitos de condições de contorno
        self._check_cc()

    @property
    def carga(self):
        """Cargas nodais externas por componente global.

        Exemplo: {'fx': 10.0, 'my': -5.0}
        """
        return dict(self.__carga)

    @carga.setter
    def carga(self, valores):
        """Define carregamento nodal por dicionário ou vetor.

        Exemplo: {'fx': 10.0, 'my': -5.0} ou [10.0, 0.0, -5.0]; use None para limpar.
        """
        if valores is None:
            self.__carga = {}
            return

        if isinstance(valores, dict):
            nova_carga = dict(valores)
        elif isinstance(valores, (list, tuple, np.ndarray)):
            if len(valores) != len(self.forcas_globais):
                raise ValueError("O vetor de carga tem comprimento incorreto")
            if not all(self._eh_numero(valor) for valor in valores):
                raise TypeError("Os valores da carga devem ser numéricos")
            nova_carga = dict(zip(self.forcas_globais, valores))
        else:
            raise TypeError(
                "carga deve ser None, um dicionário ou um vetor com valores numéricos"
            )

        for forca, valor in nova_carga.items():
            self.ifl(forca)
            if not self._eh_numero(valor):
                raise TypeError(
                    f"O valor da carga '{forca}' deve ser numérico"
                )

        self.__carga = nova_carga

    @property
    def p(self):
        """Vetor de cargas nodais na ordem padrão de `forcas_globais`."""
        return np.array([
            self.__carga.get(forca, 0.0)
            for forca in self.forcas_globais
        ])

    def _check_cc(self):
        """Valida conflitos entre condições de contorno no nó.

        Um mesmo gdl não pode ser simultaneamente nulo e prescrito,
        nem nulo e elástico.
        """
        desloc_nulos = set(self.deslocamentos_nulos)
        desloc_prescritos = set(self.deslocamentos_prescritos.keys())
        desloc_elasticos = set(self.apoio_elastico.keys())

        check = desloc_nulos.isdisjoint(desloc_prescritos) and desloc_nulos.isdisjoint(
            desloc_elasticos
        )
        if not check:
            raise Exception(
                "Há mais de um c.c. definida para os deslocamentos no nó"
            )

    def __str__(self):
        """Resumo textual do nó para inspeção rápida em debug/print."""
        linhas = []

        coords = ", ".join(
            [f"{eixo}: {valor}" for eixo, valor in zip(self.eixos_globais, self.coor)]
        )
        linhas.append(f"{coords} | igdl: {self.igdl}")

        if self.carga:
            linhas.append(f"carga: {self.carga}")

        if self.deslocamentos_nulos:
            linhas.append(f"deslocamentos_nulos: {self.deslocamentos_nulos}")

        if self.deslocamentos_prescritos:
            linhas.append(
                f"deslocamentos_prescritos: {self.deslocamentos_prescritos}"
            )

        if self.apoio_elastico:
            linhas.append(f"deslocamento_elastico: {self.apoio_elastico}")

        return "\n".join(linhas)

    def __repr__(self):
        """Representação curta para recriar o nó com as coordenadas."""
        coor = [float(c) for c in self.coor]
        return f"{self.__class__.__name__}({coor})"


class NoPE(NoBase):
    """Nó para pórtico espacial."""

    tipo = "PE"
    ndim = 3
    ngdl = 6
    eixos_locais = (1, 2, 3)
    eixos_globais = ("x", "y", "z")
    gdls_globais = ("ux", "uy", "uz", "rx", "ry", "rz")
    gdls_locais = ("u1", "u2", "u3", "r1", "r2", "r3")
    forcas_globais = ("fx", "fy", "fz", "mx", "my", "mz")


class NoPP(NoBase):
    """Nó para pórtico plano."""

    tipo = "PP"
    ndim = 2
    ngdl = 3
    eixos_locais = (1, 2)
    eixos_globais = ("x", "z")
    gdls_globais = ("ux", "uz", "ry")
    gdls_locais = ("u1", "u2", "r3")
    forcas_globais = ("fx", "fz", "my")


class NoTE(NoBase):
    """Nó para treliça espacial."""

    tipo = "TE"
    ndim = 3
    ngdl = 3
    eixos_locais = (1, 2, 3)
    eixos_globais = ("x", "y", "z")
    gdls_globais = ("ux", "uy", "uz")
    gdls_locais = ("u1", "u2", "u3")
    forcas_globais = ("fx", "fy", "fz")


class NoTP(NoBase):
    """Nó para treliça plana."""

    tipo = "TP"
    ndim = 2
    ngdl = 2
    eixos_locais = (1, 2)
    eixos_globais = ("x", "z")
    gdls_globais = ("ux", "uz")
    gdls_locais = ("u1", "u2")
    forcas_globais = ("fx", "fz")


class NoGR(NoBase):
    """Nó para grelha."""

    tipo = "GR"
    ndim = 2
    ngdl = 3
    eixos_locais = (1, 3)
    eixos_globais = ("x", "y")
    gdls_globais = ("uz", "rx", "ry")
    gdls_locais = ("u2", "r1", "r3")
    forcas_globais = ("fz", "mx", "my")


