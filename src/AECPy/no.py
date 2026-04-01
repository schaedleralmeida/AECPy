"""
Módulo para definição dos nós na análise estrutural
pelo Método da Rigidez Direta (MRD) no AECPy
"""

import numpy as np


class NoBase:
    """
    Classe base para nós em análise estrutural pelo Método da Rigidez Direta no AECPy

    Atributes
    ---------
    coor: numpy.ndarray
        Coordenadas que definem a posição do nó no modelo estrutural

    ngdl: int
        Número de graus de liberdade (gdl) do nó
    igdl: tuple (int)
        Índices globais dos gdl do nó
    gdlr: tuble (bool)
        Indicador dos gdl do nó com restrições ao deslocamento por apoios
    x, y, z
    """

    _tipo_numero = (int, float, np.int32, np.float64)
    _tipo_numero_inteiro = (int, np.int32)

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
                "A classe NoBase é base. Use uma subclasse como No_PE, No_PP, No_TE, No_TP ou No_GR."
            )
        if len(coor) != self.ndim:
            raise ValueError("Número incorreto de coordenadas")
        if not all([isinstance(c, self._tipo_numero) for c in coor]):
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
        """Coordenada x global do nó no modelo estrutural"""
        return self.coor[0]

    @property
    def y(self):
        """Coordenada y global do nó no modelo estrutural"""
        return self.coor[self.eixos_globais.index("y")]

    @property
    def z(self):
        """Coordenada z global do nó no modelo estrutural"""
        return self.coor[self.eixos_globais.index("z")]

    @property
    def igdl(self):
        """Índice global dos graus de liberdade do nó"""
        return list(self._igdl)

    @igdl.setter
    def igdl(self, ig):
        if len(ig) != self.ngdl:
            raise ValueError("O número de gdl está incorreto")
        if not all(
            [
                (isinstance(i, self._tipo_numero_inteiro) and i >= 0)
                for i in ig
            ]
        ):
            raise ValueError(
                f"Os índices globais dos gdl do nó devem ser números naturais: {ig}"
            )
        self._igdl = tuple(ig)

    def idl(self, deslocamento):
        """Retorna o índice local do deslocamento no nó."""
        if deslocamento not in self.gdls_globais:
            raise ValueError(f"Deslocamento inválido: {deslocamento}")
        return self.gdls_globais.index(deslocamento)

    def idg(self, deslocamento):
        """Retorna o índice global do deslocamento no nó."""
        return self._igdl[self.idl(deslocamento)]

    def ifl(self, forca):
        """Retorna o índice local da força no nó."""
        if forca not in self.forcas_globais:
            raise ValueError(f"Força inválida: {forca}")
        return self.forcas_globais.index(forca)

    def ifg(self, forca):
        """Retorna o índice global da força no nó."""
        return self._igdl[self.ifl(forca)]

    @property
    def deslocamentos_nulos(self):
        """Nomes dos gdl com apoios (deslocamentos nulos)."""
        return list(self.__deslocamentos_nulos)

    @deslocamentos_nulos.setter
    def deslocamentos_nulos(self, valores):
        """Define quais gdl têm deslocamentos nulos (apoios).
        
        Parâmetros
        ----------
        valores : list, tuple ou similar
            Nomes dos gdl com apoios. Deve conter apenas valores em self.gdls_globais.
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
        """Dicionário com os deslocamentos prescritos no nó."""
        return dict(self.__deslocamentos_prescritos)

    @deslocamentos_prescritos.setter
    def deslocamentos_prescritos(self, valores):
        """Define os deslocamentos prescritos no nó.

        Parâmetros
        ----------
        valores : dict
            Dicionário em que cada chave é um gdl em self.gdls_globais e cada
            valor é um número correspondente ao deslocamento prescrito.
        """
        if valores is None:
            self.__deslocamentos_prescritos = {}
            self._check_cc()
            return

        if not isinstance(valores, dict):
            raise TypeError("deslocamentos_prescritos deve ser um dicionário")

        for gdl, valor in valores.items():
            self.idl(gdl)
            if not isinstance(valor, self._tipo_numero):
                raise TypeError(
                    f"O valor do deslocamento prescrito em '{gdl}' deve ser numérico"
                )

        self.__deslocamentos_prescritos = dict(valores)

        # Verificar conflitos de condições de contorno
        self._check_cc()

    @property
    def apoio_elastico(self):
        """Dicionário com os apoios elásticos no nó."""
        return dict(self.__apoio_elastico)

    @apoio_elastico.setter
    def apoio_elastico(self, valores):
        """Define os apoios elásticos no nó.

        Parâmetros
        ----------
        valores : dict
            Dicionário em que cada chave é um gdl em self.gdls_globais e cada
            valor é um número correspondente à rigidez do apoio elástico.
        """
        if valores is None:
            self.__apoio_elastico = {}
            self._check_cc()
            return

        if not isinstance(valores, dict):
            raise TypeError("apoio_elastico deve ser um dicionário")

        for gdl, valor in valores.items():
            self.idl(gdl)
            if not isinstance(valor, self._tipo_numero):
                raise TypeError(
                    f"O valor do apoio elástico em '{gdl}' deve ser numérico"
                )

        self.__apoio_elastico = dict(valores)

        # Verificar conflitos de condições de contorno
        self._check_cc()

    @property
    def carga(self):
        """Dicionário com as cargas nodais externas."""
        return dict(self.__carga)

    @carga.setter
    def carga(self, valores):
        """Define as cargas nodais externas.

        Parâmetros
        ----------
        valores : dict, list, tuple, numpy.ndarray ou None
            `None` limpa as cargas.
            Um dicionário deve usar chaves em self.forcas_globais.
            Um vetor com comprimento igual a len(self.forcas_globais) é convertido
            para um dicionário seguindo a ordem de self.forcas_globais.
        """
        if valores is None:
            self.__carga = {}
            return

        if isinstance(valores, dict):
            nova_carga = dict(valores)
        elif isinstance(valores, (list, tuple, np.ndarray)):
            if len(valores) != len(self.forcas_globais):
                raise ValueError("O vetor de carga tem comprimento incorreto")
            if not all(isinstance(valor, self._tipo_numero) for valor in valores):
                raise TypeError("Os valores da carga devem ser numéricos")
            nova_carga = dict(zip(self.forcas_globais, valores))
        else:
            raise TypeError(
                "carga deve ser None, um dicionário ou um vetor com valores numéricos"
            )

        for forca, valor in nova_carga.items():
            self.ifl(forca)
            if not isinstance(valor, self._tipo_numero):
                raise TypeError(
                    f"O valor da carga '{forca}' deve ser numérico"
                )

        self.__carga = nova_carga

    @property
    def p(self):
        """Vetor de cargas nodais externas na ordem de forcas_globais."""
        return np.array([
            self.__carga.get(forca, 0.0)
            for forca in self.forcas_globais
        ])

    def _check_cc(self):
        """
        Levanta exceção se há mais de uma condição de contorno do tipo
        "restrito", "prescrito" ou "apoio elástico", definida simultaneamente
        para um deslocamento no nó
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
        coor = [float(c) for c in self.coor]
        return f"{self.__class__.__name__}({coor})"


class No_PE(NoBase):
    """Nó para pórtico espacial."""

    tipo = "PE"
    ndim = 3
    ngdl = 6
    eixos_locais = (1, 2, 3)
    eixos_globais = ("x", "y", "z")
    gdls_globais = ("ux", "uy", "uz", "rx", "ry", "rz")
    gdls_locais = ("u1", "u2", "u3", "r1", "r2", "r3")
    forcas_globais = ("fx", "fy", "fz", "mx", "my", "mz")


class No_PP(NoBase):
    """Nó para pórtico plano."""

    tipo = "PP"
    ndim = 2
    ngdl = 3
    eixos_locais = (1, 2)
    eixos_globais = ("x", "z")
    gdls_globais = ("ux", "uz", "ry")
    gdls_locais = ("u1", "u2", "r3")
    forcas_globais = ("fx", "fz", "my")


class No_TE(NoBase):
    """Nó para treliça espacial."""

    tipo = "TE"
    ndim = 3
    ngdl = 3
    eixos_locais = (1, 2, 3)
    eixos_globais = ("x", "y", "z")
    gdls_globais = ("ux", "uy", "uz")
    gdls_locais = ("u1", "u2", "u3")
    forcas_globais = ("fx", "fy", "fz")


class No_TP(NoBase):
    """Nó para treliça plana."""

    tipo = "TP"
    ndim = 2
    ngdl = 2
    eixos_locais = (1, 2)
    eixos_globais = ("x", "z")
    gdls_globais = ("ux", "uz")
    gdls_locais = ("u1", "u2")
    forcas_globais = ("fx", "fz")


class No_GR(NoBase):
    """Nó para grelha."""

    tipo = "GR"
    ndim = 2
    ngdl = 3
    eixos_locais = (1, 3)
    eixos_globais = ("x", "y")
    gdls_globais = ("uz", "rx", "ry")
    gdls_locais = ("u2", "r1", "r3")
    forcas_globais = ("fz", "mx", "my")


# Alias para compatibilidade com código educacional
No = NoBase
