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


def identify_t3ss_genes(degs_df):
    """
    Identifica candidatos relacionados con el sistema de
    secreción tipo III mediante sus anotaciones funcionales.
    """

    # Patrones utilizados en la búsqueda original
    t3ss_pattern = (
        r"type III secretion|"
        r"type III effector|"
        r"\bHrp[A-Za-z0-9]*\b|"
        r"\bHrc[A-Za-z0-9]*\b|"
        r"\bHop[A-Za-z0-9]*\b|"
        r"\bAvr[A-Za-z0-9]*\b"
    )

    t3ss_mask = degs_df["final_product"].str.contains(
        t3ss_pattern,
        case=False,
        na=False,
        regex=True
    )

    return degs_df.loc[t3ss_mask].copy()


def classify_t3ss_genes(t3ss_genes_df):
    """
    Clasifica los candidatos T3SS según su anotación funcional.

    Los genes sin categoría o con categorías ambiguas se
    mantienen identificados para su revisión manual.
    """

    # Patrones de clasificación del análisis original
    category_patterns = {
        "Efectores": (
            r"\beffector\b|"
            r"\bHop[A-Za-z0-9-]*\b|"
            r"\bAvr[A-Za-z0-9-]*\b|"
            r"\bXop[A-Za-z0-9-]*\b"
        ),
        "Componentes estructurales y exportacion": (
            r"type III secretion (?:system )?"
            r"(?:protein|subunit|ATPase|helper|major pilin)|"
            r"type III secretion system "
            r"(?:inner|outer|export|stator|gatekeeper)|"
            r"\bHrp pili protein\b|"
            r"\bHrpF\b|"
            r"\bHrpT\b"
        ),
        "Chaperonas": (
            r"type III secretion system chaperone|"
            r"CesT family"
        ),
        "Reguladores": (
            r"\bHrpL\b|"
            r"\bHrpR\b|"
            r"\bHrpS\b|"
            r"\bHrpV\b"
        )
    }

    def assign_category(product):
        # Identificar las categorías que coinciden
        matches = [
            category
            for category, pattern in category_patterns.items()
            if re.search(pattern, str(product), re.IGNORECASE)
        ]

        # Evitar asignaciones ambiguas automáticamente
        if len(matches) == 1:
            return matches[0]

        return "Revisar manualmente"

    classified_df = t3ss_genes_df.copy()

    classified_df["t3ss_category"] = (
        classified_df["final_product"].apply(assign_category)
    )

    return classified_df


def summarize_t3ss_genes(classified_df):
    """
    Resume la cantidad de genes inducidos y reprimidos
    por categoría funcional del T3SS.
    """

    summary_df = pd.crosstab(
        classified_df["t3ss_category"],
        classified_df["regulacion"]
    )

    return summary_df

