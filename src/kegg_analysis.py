
"""
Funciones para consultar y procesar asociaciones gen-ruta de KEGG.

Organismo de estudio:
Pseudomonas syringae pv. tomato DC3000.
"""

from io import StringIO
from urllib.error import HTTPError, URLError
from Bio.KEGG import REST
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import pandas as pd


def fetch_kegg_associations(organism_code="pst"):
    """
    Consulta KEGG para obtener asociaciones entre genes y rutas del organismo.
    Utiliza el servicio KEGG REST y transforma la respuesta tabulada
    en una tabla con identificadores de genes y rutas metabólicas.

    Parameters
    ----------
    organism_code : str, optional
        Código de organismo reconocido por KEGG. Por defecto, 'pst',
        correspondiente a Pseudomonas syringae pv. tomato DC3000.

    Returns
    -------
    pandas.DataFrame
        Asociaciones con las columnas 'kegg_gene' y 'kegg_pathway'.

    Raises
    ------
    TypeError
        Si organism_code no es una cadena de texto.
    ValueError
        Si el código está vacío o KEGG no devuelve asociaciones.
    ConnectionError
        Si la consulta a KEGG falla por un error HTTP o de conexión.
    """

    if not isinstance(organism_code, str):
        raise TypeError("organism_code debe ser una cadena de texto.")
    if not organism_code.strip():
        raise ValueError("organism_code no puede estar vacío.")

    try:
        response = REST.kegg_link("pathway", organism_code)
        content = response.read()
    except (HTTPError, URLError) as exc:
        raise ConnectionError(
            f"No fue posible consultar las asociaciones KEGG de "
            f"{organism_code}. Verifica la conexión y el código del organismo."
        ) from exc

    if not content.strip():
        raise ValueError(
            f"KEGG no devolvió asociaciones para {organism_code}."
        )
    associations_df = pd.read_csv(
        StringIO(content),
        sep="\t",
        header=None,
        names=["kegg_gene", "kegg_pathway"]
    )

    return associations_df


def normalize_kegg_associations(associations_df):
    """
    Normaliza los identificadores de genes y rutas obtenidos de KEGG.
    Por ejemplo, 'pst:PSPTO_0083' se convierte en 'PSPTO0083' y
    'path:pst03010' en 'pst03010'. Conserva las demás columnas.

    Parameters
    ----------
    associations_df : pandas.DataFrame
        Tabla de asociaciones con las columnas 'kegg_gene' y
        'kegg_pathway', como la obtenida de fetch_kegg_associations().

    Returns
    -------
    pandas.DataFrame
        Copia de la tabla con la columna 'locus_key' y los
        identificadores de rutas normalizados.

    Raises
    ------
    TypeError
        Si associations_df no es un DataFrame o sus columnas de
        identificadores no contienen cadenas de texto o valores nulos.
    ValueError
        Si faltan las columnas necesarias para la normalización.
    """

    if not isinstance(associations_df, pd.DataFrame):
        raise TypeError("associations_df debe ser un DataFrame.")

    required = {"kegg_gene", "kegg_pathway"}
    missing = required - set(associations_df.columns)
    if missing:
        raise ValueError(
            f"Faltan columnas en associations_df: {sorted(missing)}. "
            "Solución: utiliza la salida de fetch_kegg_associations()."
        )

    for column in required:
        if not associations_df[column].dropna().map(
            lambda value: isinstance(value, str)
        ).all():
            raise TypeError(
                f"La columna '{column}' debe contener identificadores de texto."
            )

    normalized_df = associations_df.copy()

    normalized_df["locus_key"] = (
        normalized_df["kegg_gene"]
        .str.replace(r"^[^:]+:", "", regex=True)
        .str.replace("_", "", regex=False)
    )
    normalized_df["kegg_pathway"] = (
        normalized_df["kegg_pathway"]
        .str.replace(r"^path:", "", regex=True)
    )

    return normalized_df


def fetch_kegg_pathway_names(organism_code="pst"):
    """
    Obtiene el catálogo de rutas KEGG y sus nombres para un organismo.
    Consulta KEGG REST y prepara una tabla con los identificadores
    normalizados de las rutas y sus descripciones funcionales.

    Parameters
    ----------
    organism_code : str, optional
        Código de organismo de KEGG. Por defecto, 'pst', para
        Pseudomonas syringae pv. tomato DC3000.

    Returns
    -------
    pandas.DataFrame
        Tabla con las columnas 'kegg_pathway' y 'pathway_name'.

    Raises
    ------
    TypeError
        Si organism_code no es una cadena de texto.
    ValueError
        Si el código está vacío o KEGG no devuelve rutas.
    ConnectionError
        Si la consulta a KEGG falla por un error HTTP o de conexión.
    """

    if not isinstance(organism_code, str):
        raise TypeError("organism_code debe ser una cadena de texto.")
    if not organism_code.strip():
        raise ValueError("organism_code no puede estar vacío.")

    # Consultar el catálogo de rutas del organismo
    try:
        response = REST.kegg_list("pathway", organism_code)
        content = response.read()
    except (HTTPError, URLError) as exc:
        raise ConnectionError(
            f"No fue posible consultar las rutas KEGG de "
            f"{organism_code}. Verifica la conexión y el código del organismo."
        ) from exc

    if not content.strip():
        raise ValueError(f"KEGG no devolvió rutas para {organism_code}.")

    # Convertir la respuesta en un DataFrame
    pathway_names_df = pd.read_csv(
        StringIO(content),
        sep="\t",
        header=None,
        names=["kegg_pathway", "pathway_name"]
    )

    # Normalizar los identificadores de las rutas
    pathway_names_df["kegg_pathway"] = (
        pathway_names_df["kegg_pathway"]
        .str.replace(r"^path:", "", regex=True)
    )

    return pathway_names_df


def run_kegg_enrichment(target_genes, pathway_universe_df):
    """
    Evalúa la sobrerrepresentación de genes de interés en rutas KEGG.
    Realiza un análisis ORA mediante la cola superior de una
    distribución hipergeométrica. Considera como universo los genes
    únicos presentes en pathway_universe_df y conserva únicamente los
    genes de interés que pertenecen a ese universo. Corrige los valores
    p de las rutas evaluadas con el método de Benjamini-Hochberg.

    Parameters
    ----------
    target_genes : iterable of str
        Identificadores de genes de interés. Se aceptan listas,
        conjuntos, índices de pandas y otros iterables de cadenas.
    pathway_universe_df : pandas.DataFrame
        Asociaciones gen-ruta del universo de referencia, con las
        columnas 'locus_key', 'kegg_pathway' y 'pathway_name'.

    Returns
    -------
    pandas.DataFrame
        Tabla ordenada por 'padj' con identificador y nombre de ruta,
        cantidad de genes de la ruta, solapamiento con genes de interés,
        valor p, genes compartidos y valor p ajustado.

    Raises
    ------
    TypeError
        Si target_genes no es un iterable de identificadores de texto,
        o pathway_universe_df no es un DataFrame.
    ValueError
        Si faltan columnas, no hay genes en el universo de referencia
        o no existen rutas KEGG evaluables.
    """

    if not isinstance(pathway_universe_df, pd.DataFrame):
        raise TypeError("pathway_universe_df debe ser un DataFrame.")

    required = {"locus_key", "kegg_pathway", "pathway_name"}
    missing = required - set(pathway_universe_df.columns)
    if missing:
        raise ValueError(
            f"Faltan columnas en pathway_universe_df: {sorted(missing)}."
        )

    if isinstance(target_genes, (str, bytes)):
        raise TypeError(
            "target_genes debe ser un iterable de identificadores, "
            "no una cadena de texto individual."
        )
    try:
        target_genes = set(target_genes)
    except (TypeError, ValueError) as exc:
        raise TypeError(
            "target_genes debe ser un iterable de identificadores válidos."
        ) from exc

    if any(not isinstance(gene, str) for gene in target_genes):
        raise TypeError("target_genes debe contener identificadores de texto.")

    for column in required:
        if not pathway_universe_df[column].dropna().map(
            lambda value: isinstance(value, str)
        ).all():
            raise TypeError(
                f"La columna '{column}' debe contener identificadores de texto."
            )

    if pathway_universe_df["locus_key"].dropna().empty:
        raise ValueError("El universo de referencia no contiene genes válidos.")

    if pathway_universe_df["kegg_pathway"].dropna().empty:
        raise ValueError("El universo de referencia no contiene rutas KEGG.")

    if pathway_universe_df.loc[
        pathway_universe_df["kegg_pathway"].notna(), "locus_key"
    ].isna().any():
        raise ValueError(
            "Hay rutas KEGG asociadas a genes sin 'locus_key'. "
            "Solución: revisa las asociaciones del universo de referencia."
        )

    # Obtener los genes unicos del universo de referencia
    background_genes = set(
        pathway_universe_df["locus_key"].dropna()
    )

    # Conservar los genes de interes presentes en el universo
    target_genes = set(target_genes) & background_genes

    population_size = len(background_genes)
    sample_size = len(target_genes)

    results = []

    # Evaluar cada ruta KEGG
    for pathway_id, group_df in pathway_universe_df.groupby(
        "kegg_pathway"
    ):

        # Identificar los genes asociados con la ruta
        pathway_genes = set(group_df["locus_key"])

        pathway_size = len(pathway_genes)

        # Identificar los genes de interes presentes en la ruta
        overlapping_genes = target_genes & pathway_genes
        overlap_size = len(overlapping_genes)

        # Probabilidad de observar al menos cantidad de genes mediante distribucion hipergeometrica
        pvalue = hypergeom.sf(
            overlap_size - 1,
            population_size,
            pathway_size,
            sample_size
        )
        results.append({
            "kegg_pathway": pathway_id,
            "pathway_name": group_df["pathway_name"].iloc[0],
            "pathway_gene_count": pathway_size,
            "overlap_count": overlap_size,
            "pvalue": pvalue,
            "overlap_genes": ", ".join(sorted(overlapping_genes))
        })

    # Construir la tabla de resultados
    results_df = pd.DataFrame(results)

    # Corregir los valores p por pruebas multiples
    results_df["padj"] = multipletests(results_df["pvalue"], method="fdr_bh")[1]

    # Ordenar las rutas por significancia estadistica
    results_df = results_df.sort_values("padj").reset_index(drop=True)

    return results_df

