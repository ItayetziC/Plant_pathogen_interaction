
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
