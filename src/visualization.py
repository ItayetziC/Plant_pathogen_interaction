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
