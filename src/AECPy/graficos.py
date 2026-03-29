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
    Dx = np.ptp(x)
    Dy = np.ptp(y)
    Dz = np.ptp(z)
    if Dx == 0: Dx = 0.1 * max(Dx, Dy, Dz)
    if Dy == 0: Dy = 0.1 * max(Dx, Dy, Dz)
    if Dz == 0: Dz = 0.1 * max(Dx, Dy, Dz)

    ax.set_box_aspect(
        (folga * Dx, folga * Dy, folga * Dz)
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

def modelo_2d(nos, els, res_els=None, tipo_diagrama=None, mostra_forcas=True):
    """
    Cria uma figura para a visualização 2D do modelo estrutural e, opcionalmente, de diagramas de esforços.

    Parâmetros
    ----------
    nos: list (of No)
        Nós do modelo
    els: list (of Elemento)
        Elementos do modelo
    res_els: list (of dict), opcional
        Resultados dos elementos (para diagramas de esforços)
    tipo_diagrama: str, opcional
        Tipo de diagrama a desenhar: 'fletor', 'cortante', 'normal' ou None
    mostra_forcas: bool, opcional
        Se True, desenha as forças aplicadas nos nós
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
    Dz = np.ptp(z)
    Dx = np.ptp(x)
    if Dz == 0: Dz =0.1*Dx
    if Dx == 0: Dx = 0.1*Dz 

    ax.set_box_aspect( Dz / Dx )

    # desenhando os elementos
    estilo_els_mod = dict(estilo_elementos)
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

    # Desenhar diagramas de esforços, se solicitado
    if tipo_diagrama is not None and res_els is not None:
        xlim = ax.get_xlim()
        zlim = ax.get_ylim()
        Dx = xlim[1] - xlim[0]
        Dz = zlim[1] - zlim[0]
        escala = 0.1 * min(Dx, Dz)
        if tipo_diagrama == 'fletor':
            chave = 'M3'
            cor = '#C54E4E'
        elif tipo_diagrama == 'cortante':
            chave = 'V2'
            cor = '#312EC9'
        elif tipo_diagrama == 'normal':
            chave = 'N'
            cor = '#55CF88'
        else:
            raise ValueError("tipo_diagrama deve ser 'fletor', 'cortante', 'normal' ou None")

        # Encontrar o maior valor absoluto para normalizar
        max_x = max(max(abs(res_el[chave])) for res_el in res_els)
        for i, el in enumerate(els):
            x_m = res_els[i]['x']
            z_m = res_els[i][chave] * escala / max_x
            x_r, z_r = rotacionar_dados_vetor(x_m, z_m, el.e1, x0=0, z0=0)
            x_rt, z_rt = transladar_dados(x_r, z_r, dx=el.noI.x, dz=el.noI.z)
            estilo_els_mod["color"] = cor
            estilo_els_mod["linewidth"] = 4
            ax.plot(x_rt, z_rt, **estilo_els_mod)

            # Linhas tracejadas conectando extremos
            x_d = [el.noI.x, x_rt[0]]
            z_d = [el.noI.z, z_rt[0]]
            ax.plot(x_d, z_d, linestyle='dashed', color=cor)
            x_d = [el.noJ.x, x_rt[-1]]
            z_d = [el.noJ.z, z_rt[-1]]
            ax.plot(x_d, z_d, linestyle='dashed', color=cor)

    # desenhando os nós
    ax.scatter(x, z, **estilo_nos)

    xlim = ax.get_xlim()
    zlim = ax.get_ylim()
    Dx = xlim[1] - xlim[0]
    Dz = zlim[1] - zlim[0]
    escala_seta = 0.1 * min(Dx, Dz)

    if mostra_forcas:
        for no in nos:
            desenhar_forcas_no(ax, no, escala_seta)

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


def unir_rel(res_els:list[dict]) -> dict:
    """
    Une os resultados de vários elementos em um único dicinário acumulando a coordenada x
    e concatenando os valores das outras quantidades
    """
    res_uni = dict(res_els[0])  # inicia com o primeiro elemento

    for res in res_els[1:]:
        for key in res_uni.keys():
           val = res[key] if key != 'x' else res[key] + res_uni['x'][-1]
           res_uni[key] = np.concatenate((res_uni[key], val))

    return res_uni

def desenhar_forcas_no(ax, no, escala_seta):
    """
    Desenha setas de forças e arcos de momento aplicados em um nó no gráfico ax.

    Parâmetros
    ----------
    ax : matplotlib.axes.Axes
        Eixo onde desenhar.
    no : objeto No
        Nó com atributos .x, .z, .p e .info()['forças globais'].
    escala_seta : float
        Escala para o tamanho das setas/arcos.
    """
    if hasattr(no, 'p') and len(no.p) > 0:
        for i, forca in enumerate(no.info()['forças globais']):
            if no.p[i] != 0:
                x0 = no.x
                z0 = no.z
                # Desenha seta para força fx ou fz
                if forca == 'fx' or forca == 'fz':
                    if forca == 'fx':
                        # Seta na direção x
                        x1 = x0 + escala_seta
                        z1 = z0
                        x2 = x1
                        z2 = z1
                    elif forca == 'fz':
                        # Seta na direção z
                        x1 = x0
                        z1 = z0 + escala_seta
                        x2 = x0 + escala_seta/2
                        z2 = z0 + escala_seta/2

                    # Desenha a seta no gráfico
                    ax.annotate('', xy=(x0, z0), xytext=(x1, z1),
                        arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=8),
                        annotation_clip=False)
                    # Adiciona o valor da força como texto próximo à seta
                    ax.annotate(f'{forca} = {-1*no.p[i]:.2f}', xy=(x0, z0), xytext=(x2, z2),
                        annotation_clip=False)

                # Desenha arco e seta para momento my
                elif forca == 'my':
                    import matplotlib.patches as mpatches
                    # Calcula pontos para desenhar o arco e a seta do momento
                    x2 = x0 + escala_seta * np.cos(np.radians(70))/2
                    z2 = z0 + escala_seta * np.sin(np.radians(70))/2
                    x3 = x0 + escala_seta * np.cos(np.radians(30))/2
                    z3 = z0 + escala_seta * np.sin(np.radians(30))/2
                    x1 = x0 + escala_seta/2
                    z1 = z0 - escala_seta/2

                    # Desenha o arco representando o momento
                    arc = mpatches.Arc((x0, z0), escala_seta, escala_seta, angle=0, theta1=215, theta2=30,
                                    color='black', linewidth=2, linestyle='-')
                    ax.add_patch(arc)

                    # Desenha a seta na extremidade do arco
                    ax.annotate('', xy=(x2, z2), xytext=(x3, z3),
                                arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=8),
                                annotation_clip=False)
                    # Adiciona o valor do momento como texto próximo ao arco
                    ax.annotate(f'{forca} = {no.p[i]:.2f}', xy=(x1, z1),
                            annotation_clip=False)

def rotacionar_dados_vetor(x, z, v_unit, x0=0, z0=0):
    """
    Rotaciona uma série de dados (x, z) em torno do ponto (x0, z0) pelo ângulo entre o vetor unitário v_unit e a horizontal (eixo x), no sentido anti-horário.

    Parâmetros
    ----------
    x : array-like
        Coordenadas x dos pontos.
    z : array-like
        Coordenadas z dos pontos.
    v_unit : array-like
        Vetor unitário (vx, vz) que define o ângulo de rotação em relação à horizontal.
    x0, z0 : float, optional
        Centro de rotação. Padrão é a origem (0,0).

    Retorna
    -------
    x_rot, z_rot : np.ndarray
        Novas coordenadas rotacionadas.
    """
    x = np.asarray(x)
    z = np.asarray(z)
    vx, vz = v_unit
    # Calcula o ângulo entre o vetor e a horizontal (eixo x)
    ang_rad = np.arctan2(vz, vx)
    # Translada para a origem
    x_t = x - x0
    z_t = z - z0
    # Matriz de rotação
    cos_a = np.cos(ang_rad)
    sin_a = np.sin(ang_rad)
    x_rot = cos_a * x_t - sin_a * z_t + x0
    z_rot = sin_a * x_t + cos_a * z_t + z0
    return x_rot, z_rot

def transladar_dados(x, z, dx=0, dz=0):
    """
    Translata uma série de dados (x, z) por (dx, dz).

    Parâmetros
    ----------
    x : array-like
        Coordenadas x dos pontos.
    z : array-like
        Coordenadas z dos pontos.
    dx : float, optional
        Deslocamento em x. Padrão é 0.
    dz : float, optional
        Deslocamento em z. Padrão é 0.

    Retorna
    -------
    x_t, z_t : np.ndarray
        Novas coordenadas transladadas.
    """
    x = np.asarray(x)
    z = np.asarray(z)
    x_t = x + dx
    z_t = z + dz
    return x_t, z_t