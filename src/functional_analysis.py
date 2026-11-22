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

def extract_genomic_region(
    genbank_file,
    start_locus=1369,
    end_locus=1406
):
    """
    Extrae los genes CDS de una región del cromosoma
    utilizando sus identificadores PSPTO.

    Se consideran tanto locus_tag como old_locus_tag
    para recuperar los identificadores históricos.
    """

    # Leer el archivo GenBank del cromosoma
    record = SeqIO.read(genbank_file, "genbank")

    genomic_genes = []

    # Recorrer las características anotadas
    for feature in record.features:

        if feature.type != "CDS":
            continue

        current_locus = feature.qualifiers.get(
            "locus_tag", [""]
        )[0]

        previous_loci = feature.qualifiers.get(
            "old_locus_tag", []
        )

        # Revisar identificadores actuales y anteriores
        locus_ids = [current_locus] + previous_loci

        for locus_id in locus_ids:

            match = re.fullmatch(
                r"PSPTO(\d+)",
                locus_id
            )

            if match is None:
                continue

            locus_number = int(match.group(1))

            if start_locus <= locus_number <= end_locus:

                genomic_genes.append({
                    "locus_key": locus_id,
                    "current_locus": current_locus,
                    "start": int(feature.location.start) + 1,
                    "end": int(feature.location.end),
                    "strand": feature.location.strand,
                    "genbank_product": feature.qualifiers.get(
                        "product", ["Sin anotación"]
                    )[0]
                })

                break

    # Construir la tabla del contexto genómico
    genomic_context_df = pd.DataFrame(
        genomic_genes,
        columns=[
            "locus_key",
            "current_locus",
            "start",
            "end",
            "strand",
            "genbank_product"
        ]
    )

    if genomic_context_df.empty:
        raise ValueError(
            "No se encontraron genes en la región solicitada."
        )

    return genomic_context_df.sort_values(
        "start"
    ).reset_index(drop=True)


def prepare_genomic_context(
    genomic_context_df,
    expression_df,
    target_locus="PSPTO1385",
    flank_genes=3
):
    """
    Integra las coordenadas genómicas de un gen y sus vecinos
    con los resultados de expresión diferencial.

    Selecciona los genes vecinos según su posición cromosómica.
    """

    # Ordenar los genes por sus coordenadas
    context_df = genomic_context_df.sort_values(
        "start"
    ).reset_index(drop=True)

    # Localizar el gen de interés
    positions = context_df.index[
        context_df["locus_key"] == target_locus
    ].tolist()

    if not positions:
        raise ValueError(
            f"No se encontró el gen {target_locus}."
        )

    target_position = positions[0]

    # Seleccionar los genes vecinos
    region_df = context_df.iloc[
        max(0, target_position - flank_genes):
        target_position + flank_genes + 1
    ].copy()

    # Incorporar la información de expresión diferencial
    expression_columns = [
        "locus_key",
        "log2FoldChange",
        "adj_pvalue",
        "regulacion"
    ]

    region_df = region_df.merge(
        expression_df[expression_columns],
        on="locus_key",
        how="left",
        validate="one_to_one"
    )

    # Identificar el gen de interés
    region_df["is_target"] = (
        region_df["locus_key"] == target_locus
    )

    return region_df.sort_values("start").reset_index(drop=True)