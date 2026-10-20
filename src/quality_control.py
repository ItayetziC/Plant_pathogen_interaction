"""Funciones para el control de calidad de las muestras."""

from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# mismos colores en todas las figuras para que cada condicion se reconozca
CONDITION_COLORS = {"plant": "tab:green", "culture": "tab:orange"}
CONDITION_LABELS = {"plant": "in planta", "culture": "King's B"}


def _check_conditions(meta: pd.DataFrame) -> None:
    """Comprueba que ``meta`` tenga ``condition`` y condiciones conocidas.

    Parameters
    ----------
    meta : pandas.DataFrame
        Tabla indexada por muestra.

    Raises
    ------
    ValueError
        Si falta ``condition`` o hay una condición sin color asignado.
    """
    if "condition" not in meta.columns:
        raise ValueError("meta debe tener la columna condition")
    unknown = set(meta["condition"]) - set(CONDITION_COLORS)
    if unknown:
        raise ValueError(
            f"Condiciones sin color asignado: {unknown}. Solución: revisa "
            "classify_sample() o añade el color en CONDITION_COLORS."
        )


def _save_figure(fig: plt.Figure, output_path: str | Path,
                 dpi: int = 150) -> None:
    """Guarda una figura con un mensaje claro si la carpeta falla.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figura a guardar.
    output_path : str or pathlib.Path
        Ruta del PNG.
    dpi : int, optional
        Resolución. Por defecto 150.

    Raises
    ------
    FileNotFoundError
        Si la carpeta no existe o no se puede escribir en ella.
    """
    try:
        fig.savefig(output_path, dpi=dpi)
    except (FileNotFoundError, PermissionError) as error:
        raise FileNotFoundError(
            f"No se pudo guardar {output_path}. Solución: crea la carpeta "
            "results/figures o revisa sus permisos."
        ) from error


def library_metrics(counts: pd.DataFrame,
                    meta: pd.DataFrame) -> pd.DataFrame:
    """Calcula la profundidad y la mediana de conteos de cada muestra.

    Se usa con la matriz SIN filtrar: el filtro quita justo los genes
    poco expresados y haría que las muestras parecieran más parejas.

    Parameters
    ----------
    counts : pandas.DataFrame
        Conteos crudos genes × muestras.
    meta : pandas.DataFrame
        Metadatos indexados por muestra, con la columna ``condition``.

    Returns
    -------
    pandas.DataFrame
        Una fila por muestra (en el orden de ``meta``) con las columnas
        ``depth_millions``, ``median_counts`` y ``condition``.

    Raises
    ------
    TypeError
        Si alguna entrada no es DataFrame o los conteos no son numéricos.
    ValueError
        Si las muestras no coinciden, la matriz está vacía, hay conteos
        negativos o condiciones desconocidas.

    Warns
    -----
    UserWarning
        Si alguna muestra tiene mediana 0, porque entonces la mediana no
        distingue esa muestra de las demás.
    """
    if not isinstance(counts, pd.DataFrame) or not isinstance(
        meta, pd.DataFrame
    ):
        raise TypeError("counts y meta deben ser DataFrame")
    _check_conditions(meta)
    if counts.empty:
        raise ValueError(
            "counts está vacío. Solución: genéralo con build_count_table()."
        )
    if set(counts.columns) != set(meta.index):
        raise ValueError(
            "Las muestras de counts y meta no coinciden. Solución: usa "
            "ambas tablas de build_count_table()."
        )
    non_numeric = [c for c in counts.columns
                   if not pd.api.types.is_numeric_dtype(counts[c])]
    if non_numeric:
        raise TypeError(f"Columnas no numéricas en counts: {non_numeric}")
    # un conteo negativo indica una matriz corrupta o ya transformada
    if (counts < 0).any().any():
        raise ValueError(
            "Hay conteos negativos. Solución: usa la matriz cruda, no "
            "una ya normalizada o en escala log."
        )
    # mismo orden que meta para que cada valor quede con su condicion
    counts = counts[meta.index]
    metrics = pd.DataFrame({
        "depth_millions": counts.sum() / 1e6,
        "median_counts": counts.median(),
        "condition": meta["condition"],
    })
    zero_median = list(metrics.index[metrics["median_counts"] == 0])
    if zero_median:
        warnings.warn(
            f"Mediana de conteos igual a 0 en {zero_median}: más de la "
            "mitad de sus genes no tienen lecturas.", UserWarning
        )
    return metrics


def plot_library_metrics(metrics: pd.DataFrame,
                         output_path: str | Path,
                         seed: int = 0) -> plt.Figure:
    """Dibuja la profundidad y la mediana de conteos por condición.

    Parameters
    ----------
    metrics : pandas.DataFrame
        Salida de :func:`library_metrics`.
    output_path : str or pathlib.Path
        Ruta del PNG a guardar.
    seed : int, optional
        Semilla del desplazamiento horizontal de los puntos, para que la
        figura salga igual en cada ejecución. Por defecto 0.

    Returns
    -------
    matplotlib.figure.Figure
        Figura generada.

    Raises
    ------
    TypeError
        Si ``metrics`` no es DataFrame.
    ValueError
        Si faltan columnas o hay condiciones desconocidas.
    FileNotFoundError
        Si no se puede guardar el PNG.
    """
    if not isinstance(metrics, pd.DataFrame):
        raise TypeError("metrics debe ser DataFrame")
    missing = {"depth_millions", "median_counts"} - set(metrics.columns)
    if missing:
        raise ValueError(
            f"A metrics le faltan {missing}. Solución: genérala con "
            "library_metrics()."
        )
    _check_conditions(metrics)
    panels = [
        ("depth_millions", "Profundidad\n(millones de lecturas)"),
        ("median_counts", "Mediana de conteos\npor gen"),
    ]
    # orden fijo: planta a la izquierda aunque meta venga en otro orden
    conditions = [c for c in CONDITION_COLORS
                  if c in set(metrics["condition"])]
    # un solo generador para que los puntos no dependan del panel
    rng = np.random.default_rng(seed)
    fig, axes = plt.subplots(1, len(panels), figsize=(8, 4.2))
    for ax, (column, title) in zip(axes, panels):
        for position, condition in enumerate(conditions):
            values = metrics.loc[metrics["condition"] == condition, column]
            # sin desplazamiento, muestras con valores parecidos se enciman
            jitter = rng.uniform(-0.07, 0.07, size=len(values))
            ax.scatter(position + jitter, values, s=60,
                       color=CONDITION_COLORS[condition], alpha=0.85,
                       edgecolor="white", linewidth=1, zorder=3)
            ax.hlines(values.mean(), position - 0.25, position + 0.25,
                      color=CONDITION_COLORS[condition], linewidth=2.5,
                      zorder=4)
        ax.set_xticks(range(len(conditions)))
        ax.set_xticklabels([
            f"{CONDITION_LABELS[c]}\n"
            f"n={int((metrics['condition'] == c).sum())}"
            for c in conditions
        ])
        ax.set_xlim(-0.55, len(conditions) - 0.45)
        ax.set_ylabel(title, fontsize=10.5)
        ax.grid(axis="y", alpha=0.25)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
    # titulo descriptivo: la conclusion se escribe tras ver la figura
    fig.suptitle(
        f"Profundidad y mediana de conteos de las {len(metrics)} muestras "
        , fontsize=12
    )
    fig.tight_layout()
    _save_figure(fig, output_path, dpi=200)
    return fig


def pca_samples(counts: pd.DataFrame,
                meta: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    """Calcula un PCA de las muestras con log2(CPM + 1).

    Parameters
    ----------
    counts : pandas.DataFrame
        Conteos genes × muestras.
    meta : pandas.DataFrame
        Metadatos indexados por muestra.

    Returns
    -------
    pcs : pandas.DataFrame
        Coordenadas de cada muestra en cada componente (``PC1``...).
    variance : numpy.ndarray
        Porcentaje de varianza explicada por componente.

    Raises
    ------
    TypeError
        Si alguna entrada no es DataFrame.
    ValueError
        Si las muestras no coinciden o hay columnas sin conteos.
    """
    if not isinstance(counts, pd.DataFrame) or not isinstance(
        meta, pd.DataFrame
    ):
        raise TypeError("counts y meta deben ser DataFrame")
    if set(counts.columns) != set(meta.index):
        raise ValueError(
            "Las muestras de counts y meta no coinciden. Solución: carga "
            "ambas tablas del mismo paso del notebook."
        )
    counts = counts[meta.index]
    library_size = counts.sum()
    # una muestra sin conteos daria division entre cero en el CPM
    if (library_size == 0).any():
        raise ValueError("Hay muestras con 0 conteos totales")
    # CPM quita el efecto de profundidad; log evita que dominen los genes
    # muy expresados
    log_cpm = np.log2(counts / library_size * 1e6 + 1)
    samples_by_genes = log_cpm.T
    # SVD sin centrar mediria la media, no la variacion entre muestras
    centered = samples_by_genes - samples_by_genes.mean(axis=0)
    left, singular, _ = np.linalg.svd(centered.values, full_matrices=False)
    variance = singular**2 / (singular**2).sum() * 100
    names = [f"PC{i + 1}" for i in range(len(singular))]
    pcs = pd.DataFrame(
        left * singular, index=samples_by_genes.index, columns=names
    )
    return pcs, variance


def plot_pca(pcs: pd.DataFrame, variance: np.ndarray, meta: pd.DataFrame,
             output_path: str | Path) -> plt.Figure:
    """Dibuja PC1 contra PC2, coloreado por condición.

    Parameters
    ----------
    pcs : pandas.DataFrame
        Salida de :func:`pca_samples`.
    variance : numpy.ndarray
        Varianza explicada, salida de :func:`pca_samples`.
    meta : pandas.DataFrame
        Metadatos con la columna ``condition``.
    output_path : str or pathlib.Path
        Ruta del PNG a guardar.

    Returns
    -------
    matplotlib.figure.Figure
        Figura generada.

    Raises
    ------
    ValueError
        Si falta ``condition`` o hay una condición sin color asignado.
    FileNotFoundError
        Si no se puede guardar el PNG.
    """
    _check_conditions(meta)
    fig, ax = plt.subplots(figsize=(7, 6))
    # un scatter por grupo para que la leyenda muestre cada condicion
    for condition, group in meta.groupby("condition"):
        ax.scatter(
            pcs.loc[group.index, "PC1"],
            pcs.loc[group.index, "PC2"],
            color=CONDITION_COLORS[condition],
            s=70,
            label=f"{condition} (n={len(group)})",
        )
    for name in pcs.index:
        # sin "Pto_" las etiquetas no se enciman
        ax.annotate(
            name.replace("Pto_", ""),
            (pcs.loc[name, "PC1"], pcs.loc[name, "PC2"]),
            fontsize=8,
            xytext=(4, 4),
            textcoords="offset points",
        )
    ax.set_xlabel(f"PC1 ({variance[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({variance[1]:.1f}%)")
    ax.legend()
    ax.set_title(f"PCA de las {len(pcs)} muestras (log2 CPM)")
    fig.tight_layout()
    _save_figure(fig, output_path)
    return fig
