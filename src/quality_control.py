"""Funciones para el control de calidad de las muestras (PCA)."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


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
    """
    colors = {"plant": "tab:green", "culture": "tab:orange"}
    if "condition" not in meta.columns:
        raise ValueError("meta debe tener la columna condition")
    unknown = set(meta["condition"]) - set(colors)
    if unknown:
        raise ValueError(f"Condiciones sin color asignado: {unknown}")
    fig, ax = plt.subplots(figsize=(7, 6))
    # un scatter por grupo para que la leyenda muestre cada condicion
    for condition, group in meta.groupby("condition"):
        ax.scatter(
            pcs.loc[group.index, "PC1"],
            pcs.loc[group.index, "PC2"],
            color=colors[condition],
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
    try:
        fig.savefig(output_path, dpi=150)
    except (FileNotFoundError, PermissionError) as error:
        raise FileNotFoundError(
            f"No se pudo guardar {output_path}. Solución: crea la carpeta "
            "results/figures o revisa sus permisos."
        ) from error
    return fig
