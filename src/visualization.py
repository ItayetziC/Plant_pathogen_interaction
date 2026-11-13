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
