"""
Funciones para el análisis funcional de genes diferencialmente expresados.
"""

import pandas as pd
import re
from Bio import SeqIO


def identify_unknown_function_genes(degs_df):
    """
    Identifica genes diferencialmente expresados cuya función
    es desconocida o está poco caracterizada.

    Parámetros
    ----------
    degs_df : pandas.DataFrame
        Tabla de genes diferencialmente expresados.

    Retorna
    -------
    pandas.DataFrame
        Genes con anotaciones de función desconocida o
        poco caracterizada.
    """

    # Definir los patrones de búsqueda en las anotaciones
    unknown_pattern = (
        r"hypothetical protein|"
        r"uncharacterized protein|"
        r"protein of unknown function|"
        r"\bDUF\d+\b"
    )

    # Identificar los genes que coinciden con los patrones
    unknown_mask = degs_df["final_product"].str.contains(
        unknown_pattern,
        case=False,
        na=False,
        regex=True
    )

    unknown_genes_df = degs_df.loc[unknown_mask].copy()

    return unknown_genes_df


def rank_unknown_function_genes(unknown_genes_df):
    """
    Ordena los genes de función desconocida según la magnitud
    de su cambio de expresión y su significancia estadística.

    Parámetros
    ----------
    unknown_genes_df : pandas.DataFrame
        Genes de función desconocida o poco caracterizada.

    Retorna
    -------
    pandas.DataFrame
        Genes ordenados por el valor absoluto de log2FoldChange.
    """

    # Crear una copia para evitar modificar los datos originales
    ranked_genes_df = unknown_genes_df.copy()

    # Calcular la magnitud del cambio de expresión
    ranked_genes_df["abs_log2FC"] = (
        ranked_genes_df["log2FoldChange"].abs()
    )

    # Ordenar por magnitud del cambio y valor p ajustado
    ranked_genes_df = ranked_genes_df.sort_values(
        ["abs_log2FC", "adj_pvalue"],
        ascending=[False, True]
    )

    return ranked_genes_df
