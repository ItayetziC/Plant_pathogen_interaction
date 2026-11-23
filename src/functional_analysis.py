"""
Funciones para el análisis funcional de genes diferencialmente expresados.
"""

import pandas as pd
import re
from Bio import SeqIO
from pathlib import Path


def identify_unknown_function_genes(degs_df):
    """
    Identifica genes diferencialmente expresados con funcion desconocida.
    Busca en la anotación `final_product` descripciones de proteínas
    hipoteticas, no caracterizadas o con dominios DUF. Conserva todas
    las columnas originales y devuelve una tabla vacia si no hay coincidencias.
    
    Parameters
    ----------
    degs_df : pandas.DataFrame
        Tabla de genes diferencialmente expresados anotados; debe
        contener la columna 'final_product'.
    
    Returns
    -------
    pandas.DataFrame
        Subconjunto de genes con función desconocida o poco caracterizada.
    
    Raises
    ------
    TypeError
        Si degs_df no es un DataFrame.
    ValueError
        Si falta la columna 'final_product'.
    """
    if not isinstance(degs_df, pd.DataFrame):
        raise TypeError("degs_df debe ser un DataFrame.")
    if "final_product" not in degs_df.columns:
        raise ValueError(
            "Falta la columna 'final_product' en degs_df. "
            "Solución: utiliza la tabla con anotación funcional."
        )
    
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
    Ordena genes de función desconocida por cambio de expresión.
    Calcula el valor absoluto de log2FoldChange en la columna
    'abs_log2FC' y ordena de mayor a menor magnitud. En caso de empate,
    prioriza los valores p ajustados más pequeños. No modifica la entrada.
    
    Parameters
    ----------
    unknown_genes_df : pandas.DataFrame
        Genes con función desconocida; requiere 'log2FoldChange' y
        'adj_pvalue' como columnas numéricas.
    
    Returns
    -------
    pandas.DataFrame
        Tabla ordenada, con la columna adicional 'abs_log2FC'.
    
    Raises
    ------
    TypeError
        Si la entrada no es un DataFrame o las columnas de análisis
        no son numéricas.
    ValueError
        Si faltan las columnas necesarias.
    """
    if not isinstance(unknown_genes_df, pd.DataFrame):
        raise TypeError("unknown_genes_df debe ser un DataFrame.")
    required = {"log2FoldChange", "adj_pvalue"}
    missing = required - set(unknown_genes_df.columns)
    if missing:
        raise ValueError(
            f"Faltan columnas en unknown_genes_df: {sorted(missing)}."
        )
    for column in required:
        if not pd.api.types.is_numeric_dtype(unknown_genes_df[column]):
            raise TypeError(f"La columna '{column}' debe ser numérica.")
        
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
    Identifica genes candidatos asociados al sistema de secrecion tipo III.
    Busca en 'final_product' términos relacionados con secrecion tipo III,
    Devuelve una tabla vacía si no encuentra coincidencias.
     
    Parameters
    ----------
    degs_df : pandas.DataFrame
        Tabla de genes anotados con la columna 'final_product'.
    
    Returns
    -------
    pandas.DataFrame
        Subconjunto de genes candidatos asociados al T3SS.
    
    Raises
    ------
    TypeError
        Si degs_df no es un DataFrame.
    ValueError
        Si falta la columna 'final_product'.
    """

    if not isinstance(degs_df, pd.DataFrame):
        raise TypeError("degs_df debe ser un DataFrame.")
    if "final_product" not in degs_df.columns:
        raise ValueError(
            "Falta la columna 'final_product' en degs_df. "
            "Solución: utiliza la tabla con anotación funcional."
        )

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
    Clasifica candidatos T3SS segun su anotacion funcional.
    Asigna las categorias 'Efectores', 'Componentes estructurales y
    exportacion', 'Chaperonas' o 'Reguladores' según patrones de texto.
    Las anotaciones sin coincidencia única quedan como 'Revisar manualmente'.
    Conserva la tabla original y agrega la columna 't3ss_category'.
    
    Parameters
    ----------
    t3ss_genes_df : pandas.DataFrame
        Genes candidatos T3SS con la columna 'final_product'.
    
    Returns
    -------
    pandas.DataFrame
        Tabla con la categoría funcional asignada a cada gen.
    
    Raises
    ------
    TypeError
        Si t3ss_genes_df no es un DataFrame.
    ValueError
        Si falta la columna 'final_product'.
    """

    if not isinstance(t3ss_genes_df, pd.DataFrame):
        raise TypeError("t3ss_genes_df debe ser un DataFrame.")
    if "final_product" not in t3ss_genes_df.columns:
        raise ValueError(
            "Falta la columna 'final_product' en t3ss_genes_df."
        )

    # Patrones de clasificación del analisis original
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
    Resume los genes T3SS por categoría funcional y regulación.
    Construye una tabla de contingencia que cuenta las observaciones
    de cada combinación entre 't3ss_category' y 'regulacion'.
    Las categorías ausentes no se incorporan automáticamente.
    
    Parameters
    ----------
    classified_df : pandas.DataFrame
        Tabla clasificada que contiene 't3ss_category' y 'regulacion'.
    
    Returns
    -------
    pandas.DataFrame
        Conteos de genes por categoría funcional (filas) y
        tipo de regulación (columnas).
    
    Raises
    ------
    TypeError
        Si classified_df no es un DataFrame.
    ValueError
        Si faltan columnas necesarias para el resumen.
    """
    if not isinstance(classified_df, pd.DataFrame):
        raise TypeError("classified_df debe ser un DataFrame.")
    required = {"t3ss_category", "regulacion"}
    missing = required - set(classified_df.columns)
    if missing:
        raise ValueError(
            f"Faltan columnas en classified_df: {sorted(missing)}."
        )
    
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
    """
    Extrae las CDS de una región identificada por loci PSPTO.
    Lee un registro GenBank y busca identificadores 'PSPTO' tanto en
    'locus_tag' como en 'old_locus_tag'. Conserva las coordenadas,
    orientación, identificador actual y anotación de producto de las
    CDS comprendidas en el rango inclusivo solicitado.
    
    Parameters
    ----------
    genbank_file : str or pathlib.Path
        Ruta de un archivo GenBank con un registro cromosómico.
    start_locus : int, optional
        Número inicial del rango PSPTO, inclusive; por defecto 1369.
    end_locus : int, optional
        Número final del rango PSPTO, inclusive; por defecto 1406.
    
    Returns
    -------
    pandas.DataFrame
        Tabla ordenada por coordenada inicial con las columnas
        'locus_key', 'current_locus', 'start', 'end', 'strand'
        y 'genbank_product'.
    
    Raises
    ------
    TypeError
        Si la ruta no es de tipo texto o Path, o si los límites no
        son enteros.
    FileNotFoundError
        Si el archivo GenBank no existe.
    ValueError
        Si la ruta está vacía, el rango es inválido o no se
        encuentran genes en la región solicitada.
    """

    if not isinstance(genbank_file, (str, Path)):
        raise TypeError("genbank_file debe ser una ruta str o Path.")
    if not str(genbank_file).strip():
        raise ValueError("genbank_file no puede ser una ruta vacía.")
    if not isinstance(start_locus, int) or isinstance(start_locus, bool):
        raise TypeError("start_locus debe ser un entero.")
    if not isinstance(end_locus, int) or isinstance(end_locus, bool):
        raise TypeError("end_locus debe ser un entero.")
    if start_locus < 0 or end_locus < start_locus:
        raise ValueError(
            "El rango de locus debe ser no negativo y start_locus "
            "no puede superar end_locus."
        )
    if not Path(genbank_file).is_file():
        raise FileNotFoundError(
            f"No se encontró el archivo GenBank: {genbank_file}"
        )

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
    Integra coordenadas genómicas con resultados de expresión diferencial.
    
    Parameters
    ----------
    genomic_context_df : pandas.DataFrame
        Tabla genómica con las columnas 'locus_key' y 'start', como
        la generada por extract_genomic_region().
    expression_df : pandas.DataFrame
        Resultados de expresión con 'locus_key', 'log2FoldChange',
        'adj_pvalue' y 'regulacion'.
    target_locus : str, optional
        Identificador PSPTO del gen de interés; por defecto 'PSPTO1385'.
    flank_genes : int, optional
        Número de genes vecinos a cada lado; por defecto 3.
    
    Returns
    -------
    pandas.DataFrame
        Región seleccionada, ordenada por posición cromosómica y
        enriquecida con información de expresión e 'is_target'.
    
    Raises
    ------
    TypeError
        Si las tablas, el identificador o flank_genes tienen tipos
        incorrectos.
    ValueError
        Si faltan columnas, flank_genes es negativo, el identificador
        es vacío, no se encuentra el gen objetivo o existen claves
        duplicadas incompatibles con el cruce uno a uno.
    """
    
    if not isinstance(genomic_context_df, pd.DataFrame):
        raise TypeError("genomic_context_df debe ser un DataFrame.")
    if not isinstance(expression_df, pd.DataFrame):
        raise TypeError("expression_df debe ser un DataFrame.")
    genomic_required = {"locus_key", "start"}
    expression_required = {
        "locus_key", "log2FoldChange", "adj_pvalue", "regulacion"
    }
    missing_genomic = genomic_required - set(genomic_context_df.columns)
    missing_expression = expression_required - set(expression_df.columns)
    if missing_genomic:
        raise ValueError(
            f"Faltan columnas en genomic_context_df: {sorted(missing_genomic)}."
        )
    if missing_expression:
        raise ValueError(
            f"Faltan columnas en expression_df: {sorted(missing_expression)}."
        )
    if not isinstance(target_locus, str):
        raise TypeError("target_locus debe ser una cadena de texto.")
    if not target_locus.strip():
        raise ValueError("target_locus no puede estar vacío.")
    if not isinstance(flank_genes, int) or isinstance(flank_genes, bool):
        raise TypeError("flank_genes debe ser un entero.")
    if flank_genes < 0:
        raise ValueError("flank_genes no puede ser negativo.")
    if genomic_context_df["locus_key"].duplicated().any():
        raise ValueError("genomic_context_df contiene locus_key duplicados.")
    if expression_df["locus_key"].duplicated().any():
        raise ValueError("expression_df contiene locus_key duplicados.")

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
