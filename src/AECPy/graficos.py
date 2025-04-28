"""
Este módulo implementa funções para gerar gráficos para
visualização de resuntados na análise estrutural pelo AECPy
"""

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

# from matplotlib import cm
import numpy as np
import pandas as pd

from .unidades import Conversor

__format_diagrama = {
    "size": 10,
    "color_marker": "black",
    "marker": ".",
    "color+": "b",
    "color-": "r",
    "alpha": 0.5,
    "color_hline": "black",
    "width_hline": 1.5,
}

__format_diagrama_p = {"size": 50, "color_marker": "g", "marker": "x"}


__nomes_ordenada_modificado = {
    "u1": "$u_1$",
    "u2": "$u_2$",
    "u3": "$u_3$",
    "r1": "$\\theta_1$",
    "r2": "$\\theta_2$",
    "r3": "$\\theta_3$",
    "Mt": "$M_t$",
    "M2": "$M_2$",
    "M3": "$M_3$",
    "V2": "$V_2$",
    "V3": "$V_3$",
}


def unidade_padrao_saida(quantidade):
    """
    Retorna uma string com a unidade padrão de saída para uma quantidade
    ou None se a unidade padrão não está definida para esta quantidade
    """

    _elementos = {
        "força": ["N", "V2", "V3"]
        + [f"f{i}" for i in [1, 2, 3, "x", "y", "z"]],
        "momento": ["Mt", "M2", "M3"]
        + [f"m{i}" for i in [1, 2, 3, "x", "y", "z"]],
        "comprimento": ["x", "y", "z"],
        "translação": [f"u{i}" for i in [1, 2, 3, "x", "y", "z"]],
        "rotação": [f"r{i}" for i in [1, 2, 3, "x", "y", "z"]],
    }
    _unidade_padrao_saida = {
        "força": "kN",
        "momento": "kNm",
        "comprimento": "cm",
        "translação": "mm",
        "rotação": "rad",
    }

    q = None
    if quantidade in _unidade_padrao_saida:
        q = quantidade
    else:
        for q, elementos in _elementos.items():
            if quantidade in elementos:
                break

    return _unidade_padrao_saida[q] if q else q


def diagramas(relpts, resultados=None, conv=None, eltag=None):
    """
    Cria uma figura com os diagramas de deslocamentos e esforços ao longo
    do eixos do elemento;

    Parâmetros
    ----------
    relpts: dict
        Dicionario com o valor das quantidades calculadas em pontos no elemento, onde
        key: string com o nome da quantidade (ex: 'N' representa a força normal)
        value: lista de valores da quantidade calculadas ao longo do elmento
        O dicinário deve ter uma key `x` representando as coordenadas a partir do início do elemento
        onde as quantidade forma calculadas
    resultados: list (optional)
        Nome dos resultados que deve ser usados para criar os diagramas;.
        Se não for definido, são criados diagramas para todos os resultados em relpts
    conv: Conversor (optional)
        Conversor de unidades, usado para escrever os resultados nas unidades padrão de saída
        Quando não for definido, os diagramas são criados sem a indicação de unidades.
    """

    abscissa = "x"

    if resultados:
        if abscissa not in resultados:
            resultados = ["x"] + resultados
        res = dict()
        for quant in resultados:
            res[quant] = relpts[quant]
    else:
        res = dict(relpts)

    # retirando a abscissa dos resultados e criando um dicionário.
    absc = {abscissa: res.pop(abscissa)}

    if conv:
        # converte os resultados para a unidade padrão de saída
        absc = conv_rel(absc, conv)
        res = conv_rel(res, conv)

    # criando as múltiplos subplots para cada resultado
    n_res = len(res)  # número de quantidades'
    fig, axes_res = plt.subplots(
        figsize=[7, n_res * 1.3], nrows=n_res, sharex=True
    )

    # identificando o eixo das abscissa no último subplot
    titulo_abscissa = [k for k in absc.keys()]
    titulo_abscissa = titulo_abscissa[0]
    vals_abscissa = np.array(absc[titulo_abscissa])
    axes_res[-1].set_xlabel(titulo_abscissa)

    # Criando os diagramas para cada resultado no elemento
    for (
        rr,
        ax,
    ) in zip(res.items(), axes_res):
        ax = criar_diagrama(ax, vals_abscissa, np.array(rr[1]), [rr[0], ""])

    #axes_res[-1].set_xlabel(titulo_abscissa)
    if eltag is not None:
        fig.suptitle(f'Elemento {eltag}')
    return fig


def criar_diagrama(ax, abscissa, ordenada, titulo_ordenada=["", ""]):
    """
    Cria o diagrama abscissa x ordenada no eixos ax

    ax     = o eixo do gráfico onde os pontos serão plotados

    abscissa: valores das abscissas do diagrama
        (normalmente, a distânica do ponto a partir do nó inicial)

    ordenada: valores das ordenadas do diagrama
        (deslocamento ou esforço)

    titulo_ordenada: título do eixo das ordenadas [nome,unidade]

    """
    if len(abscissa) != len(ordenada):
        raise ValueError(
            "abscissas e ordenadas com número distindo de elementos"
        )

    # insere para os valores positivos a cor color+ e para os valores negativos, a cor color-
    ax.fill_between(
        abscissa,
        ordenada,
        0,
        color=__format_diagrama["color+"],
        alpha=__format_diagrama["alpha"],
        where=(ordenada > 0),
        interpolate=True,
    )
    ax.fill_between(
        abscissa,
        ordenada,
        0,
        color=__format_diagrama["color-"],
        alpha=__format_diagrama["alpha"],
        where=(ordenada < 0),
        interpolate=True,
    )

    # plota os pontos calculados pela segmentação das barras
    ax.scatter(
        abscissa,
        ordenada,
        s=__format_diagrama["size"],
        color=__format_diagrama["color_marker"],
        marker=__format_diagrama["marker"],
        label=titulo_ordenada[0],
    )

    if titulo_ordenada[0] in __nomes_ordenada_modificado:
        titulo_eixo = __nomes_ordenada_modificado[titulo_ordenada[0]]
    else:
        titulo_eixo = titulo_ordenada[0]

    if titulo_ordenada[1] != "":
        titulo_eixo += " (" + titulo_ordenada[1] + ")"

    ax.set_ylabel(titulo_eixo)

    ax.axhline(
        0,
        color=__format_diagrama["color_hline"],
        linewidth=__format_diagrama["width_hline"],
    )
    ax.xaxis.grid()

    return ax


def modelo_3d(nos, els):
    """
    Cria uma figura para a visualização 3D do modelo estrutural

    Parâmetros
    ----------
    nos: list (of No)
        Nós do modelo
    els: list (of Elemento)
        Elementos do modelo
    """

    if nos[0].ndim != 3:
        raise ValueError("O modelo estrutural não é 3D")

    estilo_nos = {"c": "k", "s": 50, "alpha": 1.0}
    estilo_elementos = {
        "linestyle": "-",
        "color": "b",
        "linewidth": 5,
        "alpha": 1.0,
    }

    x = []
    y = []
    z = []
    for no in nos:
        x.append(no.x)
        y.append(no.y)
        z.append(no.z)

    fig = plt.figure(figsize=[10, 10])
    ax = fig.add_subplot(projection="3d", xlabel="x", ylabel="y", zlabel="z")
    # aspect ratio is 1:1:1 in data space
    folga = 1.2
    ax.set_box_aspect(
        (folga * np.ptp(x), folga * np.ptp(y), folga * np.ptp(z))
    )

    # desenhando os elementos
    estilo_els_mod = dict(estilo_elementos)
    i = 0
    secoes = set()
    for el in els:
        secoes.add(el.sec)
    secoes = list(secoes)
    legenda_sec = []
    for i, sec in enumerate(secoes):
        legenda_sec.append((f"C{i}", sec.nome))

    for el in els:
        x_el = [el.noI.x, el.noJ.x]
        y_el = [el.noI.y, el.noJ.y]
        z_el = [el.noI.z, el.noJ.z]
        i = secoes.index(el.sec)
        estilo_els_mod["color"] = legenda_sec[i][0]
        ax.plot(x_el, y_el, z_el, **estilo_els_mod)

    # desenhando os nós
    ax.scatter(x, y, z, **estilo_nos)

    patches = [mpatches.Patch(color=s[0], label=s[1]) for s in legenda_sec]
    ax.legend(handles=patches)
    # incluindo os títulos dos eixos
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    return fig



def modelo_2d(nos, els):
    """
    Cria uma figura para a visualização 2D do modelo estrutural

    Parâmetros
    ----------
    nos: list (of No)
        Nós do modelo
    els: list (of Elemento)
        Elementos do modelo
    """

    if nos[0].ndim != 2:
        raise ValueError("O modelo estrutural não é 2D")

    estilo_nos = {"c": "k", "s": 50, "alpha": 1.0}
    estilo_elementos = {
        "linestyle": "-",
        "color": "b",
        "linewidth": 5,
        "alpha": 1.0,
    }

    x = []
    z = []
    for no in nos:
        x.append(no.x)
        z.append(no.z)

    fig = plt.figure(figsize=[10, 10])
    ax = fig.add_subplot(xlabel="x", ylabel="z")
    
    ax.set_box_aspect( np.ptp(z) / np.ptp(x) )
    
    # desenhando os elementos
    estilo_els_mod = dict(estilo_elementos)
    i = 0
    secoes = set()
    for el in els:
        secoes.add(el.sec)
    secoes = list(secoes)
    legenda_sec = []
    for i, sec in enumerate(secoes):
        legenda_sec.append((f"C{i}", sec.nome))

    for el in els:
        x_el = [el.noI.x, el.noJ.x]
        z_el = [el.noI.z, el.noJ.z]
        i = secoes.index(el.sec)
        estilo_els_mod["color"] = legenda_sec[i][0]
        ax.plot(x_el, z_el, **estilo_els_mod)

    # desenhando os nós
    ax.scatter(x, z, **estilo_nos)

    patches = [mpatches.Patch(color=s[0], label=s[1]) for s in legenda_sec]
    ax.legend(handles=patches)
    # incluindo os títulos dos eixos
    ax.set_xlabel("x")
    ax.set_ylabel("z")

    return fig

def tabela_rel(relpts:dict, resultados:list[str]=[], conv=None) -> pd.DataFrame:
    """
    Cria uma tabela com os deslocamentos e esforços ao longo
    do eixos do elemento;

    Parâmetros
    ----------
    relpts: dict
        Dicionario com o valor das quantidades calculadas em pontos no elemento, onde
        key: string com o nome da quantidade (ex: 'N' representa a força normal)
        value: lista de valores da quantidade calculadas ao longo do elmento
        O dicinário deve ter uma key `x` representando as coordenadas a partir do início do elemento
        onde as quantidade forma calculadas
    resultados: list (optional)
        Nome dos resultados que deve ser usados para criar os diagramas;.
        Se não for definido, são criados diagramas para todos os resultados em relpts
    conv: Conversor (optional)
        Conversor de unidades, usado para escrever os resultados nas unidades padrão de saída
        Quando não for definido, os diagramas são criados sem a indicação de unidades.
    """

    if len(resultados)>0:
        if "x" not in resultados:
            resultados = ["x"] + resultados
        res = dict()
        for quant in resultados:
            res[quant] = relpts[quant]
    else:
        res = dict(relpts)

    if conv:
        # converte os resultados para a unidade padrão de saída
        res = conv_rel(res, conv)

    # Cria um DataFrame com os resultados
    df = pd.DataFrame(res)

    return df.T


def conv_rel(rel:dict, conv:Conversor):
    """
    Converte os resultados no elemento para as unidades de saída padrão
    e adiciona a unidade ao nome da quantidade
    """
    rel_ups = dict()
    for quant, vals in rel.items():
        qtext = quant
        ups = unidade_padrao_saida(quant)
        if ups:
            try:
                c = conv.para(ups)
                qtext += f" ({ups})"
            except Exception:
                if ups == "rad":
                    qtext += " (rad)"
                c = 1.0
            rel_ups[qtext] = [c * val for val in vals]
    return rel_ups
