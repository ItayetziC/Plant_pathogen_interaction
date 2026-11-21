
"""
Funciones para consultar y procesar asociaciones gen-ruta de KEGG.

Organismo de estudio:
Pseudomonas syringae pv. tomato DC3000.
"""

from io import StringIO
from Bio.KEGG import REST
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import pandas as pd


def fetch_kegg_associations(organism_code="pst"):
    """
    Consulta las asociaciones entre genes y rutas KEGG.
    
    Parametros
    ----------
    organismo : str
        Código KEGG del organismo.
    
    Retorna
    -------
    pandas.DataFrame
        Tabla con las columnas kegg_gene y kegg_pathway.
    """
    response = REST.kegg_link("pathway", organism_code)
    content = response.read()

    if not content.strip():
        raise ValueError(
            f"KEGG no devolvio asociaciones para {organism_code}."
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
    Normaliza los identificadores de genes y rutas KEGG.
    
    Ejemplos
    --------
    pst:PSPTO_0083 -> PSPTO0083
    path:pst03010 -> pst03010
    """
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
    Obtiene los identificadores y nombres de las rutas KEGG
    correspondientes al organismo de interés.

    Parametros
    ----------
    organism_code : str
        Codigo KEGG del organismo.

    Retorna
    -------
    pandas.DataFrame
        Tabla con los identificadores y nombres de las rutas.
    """
    # Consultar el catalogo de rutas del organismo
    response = REST.kegg_list("pathway", organism_code)
    content = response.read()

    if not content.strip():
        raise ValueError(f"KEGG no devolvio rutas para {organism_code}.")

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
    Realiza un analisis de sobrerrepresentacion (ORA) de rutas KEGG
    mediante una prueba hipergeometrica unilateral.

    Parametros
    ----------
    target_genes : iterable
        Identificadores de los genes de interes.

    pathway_universe_df : pandas.DataFrame
        Asociaciones gen-ruta de los genes evaluados en InMoose,
        incluyendo los nombres de las rutas.

    Retorna
    -------
    pandas.DataFrame
        Resultados del enriquecimiento, con valores p ajustados
        mediante el metodo de Benjamini-Hochberg.
    """

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