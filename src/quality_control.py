"""Funciones para el control de calidad de las muestras (PCA)."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def pca_samples(counts, meta):
    """Calcula un PCA de las muestras con log2(CPM + 1)."""
    # las muestras deben ser las mismas en ambas tablas
    if set(counts.columns) != set(meta.index):
        raise ValueError("Las muestras de 'counts' y 'meta' no coinciden")
    # mismo orden de columnas que los metadatos
    counts = counts[meta.index]
    # CPM corrige la profundidad; log2(CPM + 1) estabiliza la varianza
    cpm = counts / counts.sum() * 1e6
    log_cpm = np.log2(cpm + 1)
    # filas = muestras, columnas = genes
    samples_by_genes = log_cpm.T
    # resta el promedio de cada gen
    centered = samples_by_genes - samples_by_genes.mean(axis=0)
    left, singular, _ = np.linalg.svd(centered.values, full_matrices=False)
    # porcentaje de varianza que explica cada componente
    variance = singular**2 / (singular**2).sum() * 100
    names = [f"PC{i + 1}" for i in range(len(singular))]
    pcs = pd.DataFrame(
        left * singular, index=samples_by_genes.index, columns=names
    )
    return pcs, variance


def plot_pca(pcs, variance, meta, output_path):
    """Dibuja PC1 vs PC2 coloreando por condicion y rotulando muestras."""
    colors = {"plant": "tab:green", "culture": "tab:orange"}
    fig, ax = plt.subplots(figsize=(7, 6))
    # un grupo de puntos por condicion
    for condition, group in meta.groupby("condition"):
        ax.scatter(
            pcs.loc[group.index, "PC1"],
            pcs.loc[group.index, "PC2"],
            color=colors[condition],
            s=70,
            label=f"{condition} (n={len(group)})",
        )
    # rotula cada punto sin el prefijo "Pto_"
    for name in pcs.index:
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
    # el numero de muestras sale de los datos
    ax.set_title(f"PCA de las {len(pcs)} muestras (log2 CPM)")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    return fig