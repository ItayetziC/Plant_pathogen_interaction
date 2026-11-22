"""
Figuras de resultados del análisis de expresión diferencial.

Asunciones: la tabla de entrada es la que produce el notebook 02, con los
genes como índice y las columnas log2FoldChange, adj_pvalue y regulacion
(inducido, reprimido o sin cambio), clasificadas con adj_pvalue < 0.05 y
|log2FC| > 1, los mismos cortes que dibujan las líneas punteadas.

Cada función de dibujo recibe unos ejes opcionales: si no se le pasan, crea
los suyos. Así la misma función sirve sola o dentro de una figura compuesta.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from matplotlib.patches import FancyArrow, Patch

ALPHA = 0.05
LFC_CORTE = 1.0
COLUMNAS_REQUERIDAS = ["log2FoldChange", "adj_pvalue", "regulacion"]
COLORES = {
    "inducido": "#2c493b",      # verde
    "reprimido": "#591a24",     # rojo
    "sin cambio": "#8A9690",    # gris
}


def load_results(file_path: str | Path) -> pd.DataFrame:
    """
    Lee y valida la tabla de expresión diferencial.

    Parameters
    ----------
    file_path : str or pathlib.Path
        Ruta al archivo CSV con los resultados completos.

    Returns
    -------
    pandas.DataFrame
        Tabla sin los genes que quedaron con adj_pvalue faltante.

    Raises
    ------
    FileNotFoundError
        Si el archivo no existe.
    ValueError
        Si no es un CSV, si faltan columnas, si no quedan genes para
        graficar o si la tabla solo trae genes significativos.
    """
    file_path = Path(file_path)
    if file_path.suffix != ".csv":
        raise ValueError("el archivo debe ser .csv.")

    results = pd.read_csv(file_path, index_col=0)

    faltantes = set(COLUMNAS_REQUERIDAS) - set(results.columns)
    if faltantes:
        raise ValueError(f"faltan las columnas {sorted(faltantes)}.")

    # el volcano necesita la nube de genes sin cambio para tener contexto
    if (results["regulacion"] == "sin cambio").sum() == 0:
        raise ValueError(
            "la tabla no trae genes sin cambio. Solución: usa "
            "deg_results_completa.csv, no la de significativos."
        )

    # DESeq2 deja adj_pvalue en NaN para los genes que filtra; no se grafican
    results = results.dropna(subset=["adj_pvalue"])
    if results.empty:
        raise ValueError("no hay genes con adj_pvalue válido para graficar.")
    return results


def plot_volcano(results: pd.DataFrame, ax=None):
    """
    Dibuja el volcano plot sobre unos ejes.

    Parameters
    ----------
    results : pandas.DataFrame
        Tabla validada por :func:`load_results`.
    ax : matplotlib.axes.Axes, optional
        Ejes donde dibujar. Si no se pasan, se crean unos nuevos.

    Returns
    -------
    matplotlib.axes.Axes
        Los ejes con la gráfica, para seguir ajustándolos si hace falta.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 6))

    # clip evita -log10(0) = infinito en los genes más significativos
    significancia = -np.log10(results["adj_pvalue"].clip(lower=1e-300))

    sns.scatterplot(x=results["log2FoldChange"], y=significancia,
                    hue=results["regulacion"], palette=COLORES,
                    s=16, linewidth=0, alpha=0.7, ax=ax)

    ax.axhline(-np.log10(ALPHA), color="black", linestyle="--", linewidth=1)
    ax.axvline(LFC_CORTE, color="black", linestyle="--", linewidth=1)
    ax.axvline(-LFC_CORTE, color="black", linestyle="--", linewidth=1)

    ax.set_title("Expresión diferencial de DC3000 durante la infección")
    ax.set_xlabel("log2 (in planta / King's B)")
    ax.set_ylabel("-log10 (FDR)")
    ax.legend(title="")
    return ax

def plot_heatmap(log_counts: pd.DataFrame,
                 genes: pd.Index,
                 meta: pd.DataFrame,
                 ax=None,
                 etiquetas: list[str] | None = None):
    """
    Dibuja un heatmap de los genes seleccionados.

    Cada fila se estandariza por separado (puntaje z), de modo que el color
    indica cuánto se aleja cada muestra del promedio de ese gen. Sin
    estandarizar, los genes muy expresados dominarían toda la escala.

    Las columnas se ordenan por condición y las filas por log2FoldChange,
    sin agrupamiento, para que los bloques se lean de inmediato.

    Parameters
    ----------
    log_counts : pandas.DataFrame
        Expresión en log2 de conteos normalizados, genes × muestras.
    genes : pandas.Index
        Genes a graficar, en el orden deseado.
    meta : pandas.DataFrame
        Metadatos con la columna ``condition``.
    ax : matplotlib.axes.Axes, optional
        Ejes donde dibujar. Si no se pasan, se crean unos nuevos.
    etiquetas : list of str, optional
        Nombres a mostrar en las filas. Por omisión, los identificadores.

    Returns
    -------
    matplotlib.axes.Axes

    Raises
    ------
    ValueError
        Si algún gen solicitado no está en la matriz de expresión.
    """
    ausentes = [g for g in genes if g not in log_counts.index]
    if ausentes:
        raise ValueError(
            f"{len(ausentes)} genes no están en la matriz de expresión "
            f"(por ejemplo {ausentes[:3]})."
        )

    if ax is None:
        _, ax = plt.subplots(figsize=(8, 10))

    orden_muestras = meta.sort_values("condition").index
    datos = log_counts.loc[genes, orden_muestras]

    # puntaje z por fila: (valor - media del gen) / desviación del gen
    z = datos.sub(datos.mean(axis=1), axis=0).div(datos.std(axis=1), axis=0)

    sns.heatmap(z, cmap="RdBu_r", center=0, ax=ax,
                yticklabels=etiquetas if etiquetas is not None else list(genes),
                xticklabels=[s.replace("Pto_", "") for s in orden_muestras],
                cbar_kws={"label": "puntaje z", "shrink": 0.5})

    # línea que separa las dos condiciones
    n_cultivo = int((meta["condition"] == "culture").sum())
    ax.axvline(n_cultivo, color="black", linewidth=2)

    ax.set_title(f"Los {len(genes)} genes con mayor cambio de expresión")
    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.setp(ax.get_yticklabels(), fontsize=8, fontstyle="italic", rotation=0)
    plt.setp(ax.get_xticklabels(), fontsize=8, rotation=90)
    return ax


def select_top_degs(degs: pd.DataFrame, n_each: int = 20) -> pd.Index:
    """
    Elige los genes más inducidos y más reprimidos.

    Parameters
    ----------
    degs : pandas.DataFrame
        Subconjunto de genes diferencialmente expresados.
    n_each : int
        Cuántos tomar de cada dirección.

    Returns
    -------
    pandas.Index
        Identificadores ordenados de mayor a menor log2FoldChange.

    Raises
    ------
    ValueError
        Si no hay suficientes genes en alguna dirección.
    """
    inducidos = degs[degs["log2FoldChange"] > 0]
    reprimidos = degs[degs["log2FoldChange"] < 0]
    if len(inducidos) < n_each or len(reprimidos) < n_each:
        raise ValueError(
            f"no hay {n_each} genes en cada dirección "
            f"({len(inducidos)} inducidos, {len(reprimidos)} reprimidos)."
        )

    seleccion = pd.concat([
        inducidos.nlargest(n_each, "log2FoldChange"),
        reprimidos.nsmallest(n_each, "log2FoldChange"),
    ])
    return seleccion.sort_values("log2FoldChange", ascending=False).index

def plot_candidatos(results: pd.DataFrame,
                    ax=None,
                    n: int = 15,
                    base_mean_minimo: float = 10.0,
                    solo_inducidos: bool = True,
                    patron: str = "hypothetical"):
    """
    Grafica los genes sin función asignada que cambian durante la infección.
 
    Responde la pregunta 3: entre los genes anotados solo como proteína
    hipotética, cuáles muestran un patrón de expresión que los vuelve
    candidatos a caracterización funcional.
 
    Se exige una expresión mínima además de la significancia, porque un
    cambio grande en un gen apenas detectable no es evidencia sólida.
 
    Parameters
    ----------
    results : pandas.DataFrame
        Tabla anotada, con ``final_product``, ``log2FoldChange``,
        ``baseMean`` y ``regulacion``.
    ax : matplotlib.axes.Axes, optional
        Ejes donde dibujar. Si no se pasan, se crean unos nuevos.
    n : int
        Cuántos candidatos mostrar.
    base_mean_minimo : float
        Expresión media mínima para considerar al gen.
    solo_inducidos : bool
        Si es True, solo muestra los inducidos in planta.
    patron : str
        Texto a buscar en ``final_product`` para identificar genes sin
        función asignada.
 
    Returns
    -------
    matplotlib.axes.Axes
 
    Raises
    ------
    ValueError
        Si ningún gen cumple los criterios.
    """
    sin_funcion = results["final_product"].fillna("").str.contains(
        patron, case=False, regex=False
    )
    candidatos = results[
        sin_funcion
        & (results["regulacion"] != "sin cambio")
        & (results["baseMean"] >= base_mean_minimo)
    ]
    if solo_inducidos:
        candidatos = candidatos[candidatos["log2FoldChange"] > 0]
 
    if candidatos.empty:
        raise ValueError(
            f"ningún gen cumple los criterios (patrón '{patron}', "
            f"baseMean >= {base_mean_minimo})."
        )
 
    if solo_inducidos:
        candidatos = candidatos.nlargest(n, "log2FoldChange")
    else:
        # mitad y mitad, para que ninguna dirección domine la figura
        mitad = n // 2
        candidatos = pd.concat([
            candidatos[candidatos["log2FoldChange"] > 0].nlargest(mitad, "log2FoldChange"),
            candidatos[candidatos["log2FoldChange"] < 0].nsmallest(n - mitad, "log2FoldChange"),
        ])
    candidatos = candidatos.sort_values("log2FoldChange")
 
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 0.42 * len(candidatos) + 2))
 
    colores = [COLORES["inducido"] if v > 0 else COLORES["reprimido"]
               for v in candidatos["log2FoldChange"]]
    posiciones = np.arange(len(candidatos))
    ax.barh(posiciones, candidatos["log2FoldChange"], color=colores, height=0.68)
 
    # la expresión media se anota al final de cada barra
    for y, (lfc, media) in enumerate(zip(candidatos["log2FoldChange"],
                                         candidatos["baseMean"])):
        desplazamiento = 0.12 if lfc > 0 else -0.12
        ax.text(lfc + desplazamiento, y, f"{media:.0f}", va="center",
                ha="left" if lfc > 0 else "right", fontsize=8, color="#444")
 
    ax.set_yticks(posiciones)
    ax.set_yticklabels(candidatos.index, fontsize=9, fontstyle="italic")
    ax.axvline(0, color="black", linewidth=1)
    ax.set_xlabel("log2 (in planta / King's B)")
    direccion = "inducidos" if solo_inducidos else "con cambio de expresión"
    ax.set_title(
        f"Genes sin función asignada {direccion} durante la infección\n"
        f"(FDR < {ALPHA}, expresión media >= {base_mean_minimo:.0f}; "
        f"la cifra junto a cada barra es la expresión media)",
        fontsize=10.5,
    )
    margen = candidatos["log2FoldChange"].abs().max() * 0.22
    ax.set_xlim(min(0, candidatos["log2FoldChange"].min() - margen),
                max(0, candidatos["log2FoldChange"].max() + margen))
    ax.grid(axis="x", alpha=0.25)
    ax.set_axisbelow(True)
    for lado in ("top", "right"):
        ax.spines[lado].set_visible(False)
    return ax

def plot_kegg_enrichment(results_df, title, ax=None):
    """
    Genera un dot plot de rutas KEGG significativamente enriquecidas.

    Parametros
    ----------
    results_df : pandas.DataFrame
        Resultados del analisis ORA.
    title : str
        Titulo de la grafica.
    ax : matplotlib.axes.Axes, opcional
        Eje sobre el que se dibujar la grafica.

    Retorna
    -------
    matplotlib.axes.Axes
        Eje que contiene la grafica.
    """
    # Seleccionar las rutas significativamente enriquecidas
    significant_df = results_df.loc[
        results_df["padj"] < 0.05
    ].copy()

    # Simplificar los nombres de las rutas para la visualizacion
    significant_df["short_pathway_name"] = (
        significant_df["pathway_name"]
        .str.split(" - ").str[0]
    )

    if ax is None:
        _, ax = plt.subplots(figsize=(10, 5))

    if significant_df.empty:
        ax.text(
            0.5, 0.5,
            "No se identificaron rutas significativas",
            ha="center",
            va="center",
            transform=ax.transAxes
        )
        ax.set_axis_off()
        return ax

    # Calcular la significancia en escala logaritmica
    significant_df["neg_log10_padj"] = (
        -np.log10(significant_df["padj"].clip(lower=1e-300))
    )

    # Ordenar las rutas por significancia
    significant_df = significant_df.sort_values(
        "neg_log10_padj"
    )

    # Representar las rutas mediante un grafico de puntos
    scatter = ax.scatter(
        significant_df["overlap_count"],
        significant_df["short_pathway_name"],
        s=significant_df["overlap_count"] * 12,
        c=significant_df["neg_log10_padj"],
        cmap="viridis",
        alpha=0.85,
        edgecolors="black",
        linewidths=0.5
    )

    ax.set_xlabel("Numero de DEGs asociados a la ruta")
    ax.set_ylabel("Ruta KEGG")
    ax.set_title(title)

    plt.colorbar(
        scatter,
        ax=ax,
        label="-log10(FDR)"
    )

    return ax


def plot_t3ss_categories(summary_df, ax=None):
    """
    Genera una gráfica de barras horizontales apiladas
    para representar la clasificación funcional del T3SS.

    Parámetros
    ----------
    summary_df : pandas.DataFrame
        Tabla con el número de genes inducidos y reprimidos
        por categoría funcional.

    ax : matplotlib.axes.Axes, opcional
        Eje donde se dibujará la gráfica.

    Retorna
    -------
    matplotlib.axes.Axes
        Eje con la gráfica generada.
    """

    # Crear una copia de la tabla
    plot_df = summary_df.copy()

    # Asegurar que ambas categorías de regulación estén presentes
    for column in ["inducido", "reprimido"]:
        if column not in plot_df.columns:
            plot_df[column] = 0

    # Ordenar las categorías funcionales
    category_order = [
        "Efectores",
        "Componentes estructurales y exportacion",
        "Chaperonas",
        "Reguladores"
    ]

    plot_df = plot_df.reindex(category_order, fill_value=0)

    # Crear la figura cuando no se proporciona un eje
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 5))

    # Dibujar las barras horizontales apiladas
    ax.barh(
        plot_df.index,
        plot_df["inducido"],
        label="Inducidos",
        color="#4C78A8"
    )

    ax.barh(
        plot_df.index,
        plot_df["reprimido"],
        left=plot_df["inducido"],
        label="Reprimidos",
        color="#E45756"
    )

    # Agregar el total de genes al final de cada barra
    for index, (_, row) in enumerate(plot_df.iterrows()):
        total = row["inducido"] + row["reprimido"]

        ax.text(
            total + 0.3,
            index,
            str(int(total)),
            va="center"
        )

    # Configurar las etiquetas
    ax.set_xlabel("Número de genes diferencialmente expresados")
    ax.set_ylabel("Categoría funcional")
    ax.set_title("Clasificación funcional de genes relacionados con T3SS")

    ax.legend(title="Regulación")
    ax.invert_yaxis()
    ax.set_xlim(0, plot_df.sum(axis=1).max() + 4)

    return ax


def plot_genomic_context(context_df, target_locus="PSPTO1385", ax=None):
    """
    Representa el contexto genómico de un candidato mediante
    flechas que indican la posición, orientación y regulación.
    """

    if context_df.empty:
        raise ValueError("No hay genes para representar.")

    # Crear la figura
    if ax is None:
        _, ax = plt.subplots(figsize=(13, 4.5))

    # Definir los colores según la regulación
    regulation_colors = {
        "inducido": "#D95F02",
        "reprimido": "#1B9E77",
        "sin cambio": "#A0A0A0"
    }

    region_start = context_df["start"].min()
    region_end = context_df["end"].max()

    # Dibujar la línea que representa el cromosoma
    ax.plot(
        [region_start, region_end],
        [0, 0],
        color="gray",
        linewidth=1.5,
        zorder=1
    )

    # Representar cada gen
    for _, gene in context_df.iterrows():

        start = gene["start"]
        end = gene["end"]
        strand = gene["strand"]

        # Determinar la dirección de transcripción
        if strand == 1:
            x = start
            dx = end - start
        elif strand == -1:
            x = end
            dx = start - end
        else:
            continue

        regulation = gene["regulacion"]

        color = regulation_colors.get(
            regulation,
            "#A0A0A0"
        )

        is_target = gene["locus_key"] == target_locus

        # Crear la flecha que representa el gen
        arrow = FancyArrow(
            x,
            0,
            dx,
            0,
            width=0.20,
            head_width=0.32,
            head_length=min(70, abs(dx) * 0.3),
            length_includes_head=True,
            facecolor=color,
            edgecolor="black" if is_target else "none",
            linewidth=2 if is_target else 0,
            zorder=3
        )

        ax.add_patch(arrow)

        # Calcular el centro del gen
        center = (start + end) / 2

        # Mostrar el identificador
        ax.text(
            center,
            0.24,
            gene["locus_key"],
            ha="center",
            va="bottom",
            fontsize=9,
            rotation=45
        )

        # Mostrar el cambio de expresión
        if pd.notna(gene["log2FoldChange"]):
            ax.text(
                center,
                -0.28,
                f'{gene["log2FoldChange"]:+.2f}',
                ha="center",
                va="top",
                fontsize=9
            )

    # Configurar los ejes
    ax.set_xlim(region_start - 150, region_end + 150)
    ax.set_ylim(-0.65, 0.95)
    ax.set_yticks([])

    ax.set_xlabel("Posicion en el cromosoma (pb)")
    ax.set_title(
        f"Contexto genomico y expresion diferencial de {target_locus}",
        pad=18
    )

    # Crear la leyenda
    legend_items = [
        Patch(facecolor="#D95F02", label="Inducido"),
        Patch(facecolor="#1B9E77", label="Reprimido"),
        Patch(facecolor="#A0A0A0", label="Sin cambio o sin datos"),
        Patch(
            facecolor="white",
            edgecolor="black",
            linewidth=2,
            label=target_locus
        )
    ]

    ax.legend(
        handles=legend_items,
        loc="upper right",
        frameon=False
    )

    ax.ticklabel_format(style="plain", axis="x")

    return ax
