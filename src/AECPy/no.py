"""
Módulo para definição dos nós na análise estrutural
pelo Método da Rigidez Direta (MRD) no AECPy
"""

import numpy as np


class No:
    """
    Nó para análise estrutural pelo Método da Rigidez Direta no AECPy

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

    __tipo_numero = (int, float, np.int32, np.float64)
    __tipo_numero_inteiro = (int, np.int32)

    # atributos usados como default para os nós
    __igdl = ()  # índice global dos gdl do nó
    __gdlr = ()  # índice local dos graus de liberdade restritos no nó
    __gdlp = ()  # índice local dos graus de liberdade prescritos no nó
    __gdle = ()  # índice local dos graus de liberdade com apoio elástico no nó
    __carga = ()  # vetor de carga externa no nó
    __k = ()  # vetor de coeficientes de rigidez do apoio elástico no nó

    @classmethod
    def iniciar(cls, tipo) -> None:
        if tipo == "TP":  # Treliça Plana
            cls.__ndim = 2
            cls.__ngdl = 2
            cls.__eixos_locais = (1, 2)
            cls.__eixos_globais = ("x", "z")
            cls.__gdls_globais = ("ux", "uz")
            cls.__gdls_locais = ("u1", "u2")
            cls.__forcas_globais = ("fx", "fy")

        elif tipo == "TE":  # Treliça Espacial
            cls.__ndim = 3
            cls.__ngdl = 3
            cls.__eixos_locais = (1, 2, 3)
            cls.__eixos_globais = ("x", "y", "z")
            cls.__gdls_globais = ("ux", "uy", "uz")
            cls.__gdls_locais = ("u1", "u2", "u3")
            cls.__forcas_globais = ("fx", "fy", "fz")

        elif tipo == "PP":  # Pórtico Plano
            cls.__ndim = 2
            cls.__ngdl = 3
            cls.__eixos_locais = (1, 2)
            cls.__eixos_globais = ("x", "z")
            cls.__gdls_globais = ("ux", "uz", "ry")
            cls.__gdls_locais = ("u1", "u2", "r3")
            cls.__forcas_globais = ("fx", "fz", "my")

        elif tipo == "GR":  # Grelha
            cls.__ndim = 2
            cls.__ngdl = 3
            cls.__eixos_locais = (1, 3)
            cls.__eixos_globais = ("x", "y")
            cls.__gdls_globais = ("uz", "rx", "ry")
            cls.__gdls_locais = ("u2", "r1", "r3")
            cls.__forcas_globais = ("fz", "mx", "my")

        elif tipo == "PE":  # Pórtico Espacial
            cls.__ndim = 3
            cls.__ngdl = 6
            cls.__eixos_locais = (1, 2, 3)
            cls.__eixos_globais = ("x", "y", "z")
            cls.__gdls_globais = ("ux", "uy", "uz", "rx", "ry", "rz")
            cls.__gdls_locais = ("u1", "u2", "u3", "r1", "r2", "r3")
            cls.__forcas_globais = ("fx", "fy", "fz", "mx", "my", "mz")

        else:
            raise ValueError(f"O tipo `{tipo}` não é suportado")

        cls.__tipo = tipo

        cls.__info = {
            "tipo": cls.__tipo,
            "ndim": cls.__ndim,
            "ngdl": cls.__ngdl,
            "eixos locais": cls.__eixos_locais,
            "eixos globais": cls.__eixos_globais,
            "gdls globais": cls.__gdls_globais,
            "gdls locais": cls.__gdls_locais,
            "forças globais": cls.__forcas_globais,
        }

    @classmethod
    def info(cls, info_requerida=None):
        """Retorna uma informação definida na inicialização da classe"""
        if info_requerida is None:
            return dict(cls.__info)
        return cls.__info[info_requerida]

    def __init__(self, coor) -> None:
        if len(coor) != self.ndim:
            raise ValueError("Número incorreto de coordenadas")
        if not all([isinstance(c, self.__tipo_numero) for c in coor]):
            raise TypeError("As coordenadas devem ser números")
        self.__coor = np.array(coor)  # coordenadas da posição do nó

    @property
    def ndim(self):
        """Número de dimensões do espaço no modelo (2 ou 3)"""
        return self.__ndim

    @property
    def ngdl(self):
        """Número de gdl do nó"""
        return self.__ngdl

    @property
    def coor(self):
        """Coordenadas do nó"""
        return self.__coor

    @property
    def x(self):
        """Coordenada x global do nó no modelo estrutural"""
        return self.coor[0]

    @property
    def y(self):
        """Coordenada y global do nó no modelo estrutural"""
        return self.coor[self.__eixos_globais.index("y")]

    @property
    def z(self):
        """Coordenada z global do nó no modelo estrutural"""
        return self.coor[self.__eixos_globais.index("z")]

    @property
    def igdl(self):
        """Índice global dos graus de liberdade do nó"""
        return list(self.__igdl)

    @igdl.setter
    def igdl(self, ig):
        if len(ig) != self.ngdl:
            raise ValueError("O número de gdl está incorreto")
        if not all(
            [
                (isinstance(i, self.__tipo_numero_inteiro) and i >= 0)
                for i in ig
            ]
        ):
            raise ValueError(
                f"Os índices globais dos gdl do nó devem ser números naturais: {ig}"
            )
        self.__igdl = tuple(ig)

    @property
    def ngdlr(self):
        """Número de gdl restritos no nó"""
        return len(self.__gdlr)

    @property
    def il_gdlr(self):
        """Índice local dos gdl restritos, cujo deslocamento é impedido por apoios"""
        return self.__gdlr

    @property
    def igdlr(self):
        """Índice global dos gdl restritos no nó"""
        return [self.__igdl[il] for il in self.__gdlr]

    @property
    def ngdlp(self):
        """Número de gdl prescritos no nó"""
        return len(self.__gdlp)

    @property
    def il_gdlp(self):
        """Índice local dos gdl prescritos no nó"""
        return self.__gdlp

    @property
    def igdlp(self):
        """Índice global dos gdl prescritos no nó"""
        return [self.__igdl[il] for il in self.__gdlp]

    @property
    def dp(self):
        """Valor do deslocamento prescrito nos gdlp"""
        return self.__dp

    @property
    def ngdle(self):
        """Número de gdl com apoio elástico no nó"""
        return len(self.__gdle)

    @property
    def il_gdle(self):
        """Índice local dos gdl com apoio elástico no nó"""
        return self.__gdle

    @property
    def igdle(self):
        """Índice global dos gdl com apoio elástico no nó"""
        return [self.__igdl[il] for il in self.__gdle]

    @property
    def k(self):
        """Valor da constante de rigidez nos gdle"""
        return self.__k

    @property
    def p(self):
        """Vetor de força externa aplicada no nó"""
        return self.__carga

    @property
    def carregado(self):
        """Retorna True se há forças externas aplicadas no nó"""
        return len(self.__carga) > 0

    def definir_apoio(self, *args) -> None:
        """
        Define os apoios no nó, indicando os gdl com deslocamento nulo.

        Parameter
        ---------
        args: None, str
            Indicação dos gdl com deslocamento nulo no nó
            None: libera o deslocamento em todos os gdl
            'Todos': restringe o deslocamento em todos os gdl
            'ux','rz': restringe os deslocamentos nos gdl indicados e libera os demais

        Exemplo
        -------
        >>> no.definir_apoio(None)
            elimina os apoios do nó
        >>> no.definir_apoio('todos')
            restringe o deslocamento de godos os gdl do nó
        >>> no.definir_apoio( 'ux','rz')
            restringe o deslocamento dos gdl 'ux' e 'rz', liberando os demais
            e retorna um erro se 'ux' e 'rz' não são gdl do nó.
        """

        if None in args:
            self.__gdlr = ()
        elif "todos" in args:
            self.__gdlr = tuple([i for i in range(self.ngdl)])
        else:
            aux = set()
            for gdl in args:
                aux.add(self.__gdls_globais.index(gdl))
            self.__gdlr = sorted(aux)

        # verificando c.c repetidas no nó
        self._check_cc()

        return

    def definir_deslocamento(self, **kwargs):
        """
        Define os deslocamentos prescritos no nó

        Parameter
        ---------
        kwargs: dict
            key = deslocamento prescrito
            value = valor do deslocamento

        Exemplo
        -------
        >>> no.definir_deslocamento(ux=3.5, rz=1.e-3)
            define valores 3.5 e 1.e-3 para a translação na direção do eixo x (ux)
            e rotação em torno do eixo z (rz), respectivamente
        >>> no.definir_deslocamento( )
            elimina qualquer deslocamento prescrito anteriormente no nó
        """
        self.__gdlp = []
        self.__dp = []
        aux = set()
        for dir, dp in kwargs.items():
            if dir not in self.__gdls_globais:
                raise ValueError(
                    "O deslocamento {} não pode ser prescrito para {}".format(
                        dir, self.__tipo
                    )
                )
            if dir in aux:
                raise ValueError(
                    "O deslocamento {} tem mais de um valor prescritos".format(
                        dir
                    )
                )
            if not isinstance(dp, self.__tipo_numero):
                raise TypeError("O deslocamento prescrito deve ser um número")
            aux.add(dir)
            self.__gdlp.append(self.__gdls_globais.index(dir))
            self.__dp.append(dp)

        # verificando c.c repetidas no nó
        self._check_cc()

        return

    def definir_apoio_elastico(self, **kwargs):
        """
        Define os apoios elásticos no nó, indicando o gdl e a constante de rigidez do apoio
        Parameter
        ---------
        kwargs: dict
            key = deslocamento com apoio elástico
            value = constante de rigidez do apoio

        Exemplo
        -------
        >>> no.definir_apoio_elastico(ux=1.3e4, rz=2.5e5)
            define valores 1.3e4 e 2.5e5 para a constante de rigidez do apoio
            de translação na direção do eixo x (ux) e rotação em torno do eixo
            z (rz), respectivamente
        >>> no.definir_apoio_elastico( )
            elimina qualquer apoio elástico definido anteriormente no nó
        """
        self.__gdle = []
        self.__k = []
        aux = set()
        for dir, k in kwargs.items():
            if dir not in self.__gdls_globais:
                raise ValueError(
                    "O apoio elastico em {} não pode ser definido para {}".format(
                        dir, self.__tipo
                    )
                )
            if dir in aux:
                raise ValueError(
                    "O apoio elástico em {} tem mais de um valor de reigidez definido".format(
                        dir
                    )
                )
            if not isinstance(k, self.__tipo_numero):
                raise TypeError("A constante de rigidez deve ser um número")
            aux.add(dir)
            self.__gdle.append(self.__gdls_globais.index(dir))
            self.__k.append(k)

        # verificando c.c repetidas no nó
        self._check_cc()

        return

    def _check_cc(self):
        """
        Levanta exceção se há mais de uma condição de contorno do tipo
        "restrito", "prescrito" ou "apoio elástico", definida simultaneamente
        para um deslocamento no nó
        """
        # deslcamentos restritos
        rr = set(self.__gdlr)
        pp = set(self.__gdlp)
        ee = set(self.__gdle)
        check = rr.isdisjoint(pp) and rr.isdisjoint(ee) and pp.isdisjoint(ee)
        if not check:
            raise Exception(
                "Há mais de um c.c. definida para os deslocamentos no nó"
            )

    def definir_carga(self, **kwargs) -> None:
        """
        Define o carregamento no nó, indicando os gdl com força/momento externo.

        Parameter
        ---------
        kwargs: dict
            Indicação das forças e momentos aplicados no nó
            Entradas válidas:
            p=forca; define o vetor de força em todos os gdl, com len(forca)=no.ngdl.
            p=None: elimina a força em todos os gdl
            fa=forca: aplica força (float) na direção do eixo global `a`, onde `a` é x, y ou z
            ma=momento: aplica momento (float) em torno do eixo global `a`, onde `a` é x, y ou z

        Exemplo
        -------
        >>> no.definir_carga(p = None)
            elimina a carga externa no nó (no.p é [0,0,0])
        >>> no.definir_carga(p = [10,30,5])
            define uma carga externa para todos os gdl do nó, conforme ordem padrão. (no.p é [10,30,5])
        >>> no.definir_carga(fx=10, my=20)
            define força 10 no gdl correspondente a `fx` e momento 20 no gdl correspondente a `my`
            e força/momento nulos para os demais gdl. ( no.p é [10, 0, 20] se as forças permitidas são [fx,fz,my])
        """

        if "p" in kwargs:
            p = kwargs["p"]
            if p is None:
                self.__carga = ()
            elif len(p) == self.ngdl and all(
                [isinstance(v, (int, float)) for v in p]
            ):
                self.__carga = np.array(p)
            else:
                raise ValueError("O vetor de carga p está incorreto")
        else:
            self.__carga = np.zeros(self.ngdl)

            for i, p in enumerate(self.__forcas_globais):
                if p in kwargs:
                    i = self.__forcas_globais.index(p)
                    if isinstance(kwargs[p], (int, float)):
                        self.__carga[i] = kwargs[p]
                    else:
                        raise ValueError(
                            f"A carga {p} tem valor incorreto: {kwargs[p]}"
                        )

        return

    def dados(self):
        """Retorna um dicionário com dados do nó"""
        dados = dict()
        # coordandas
        for eixo, coor in zip(self.__eixos_globais, self.coor):
            dados[eixo] = coor

        # índice global do gdl
        dados["igdl"] = self.igdl

        # carregamento externo
        if self.carregado:
            for forca, val in zip(self.__forcas_globais, self.p):
                dados[forca] = val

        # Deslocamento restritos - apoios:
        if self.ngdlr > 0:
            dados["gdl_restritos"] = [
                self.__gdls_globais[i] for i in self.il_gdlr
            ]

        # Deslocamentos prescritos:
        if self.ngdlp > 0:
            for il, dp in zip(self.il_gdlp, self.dp):
                dados["dp_" + self.__gdls_globais[il]] = dp

        # Apoios elásticos:
        if self.ngdle > 0:
            for il, k in zip(self.il_gdle, self.k):
                dados["k_" + self.__gdls_globais[il]] = k

        return dados

    def __str__(self):
        dados = self.dados()
        for key, val in dados.items():
            print(f"{key}: {val}")

    def __repr__(self):
        return (
            self.__class__.__name__
            + f"( {self.coor},{self.igdl},{self.il_gdlr},{self.p})"
        )
