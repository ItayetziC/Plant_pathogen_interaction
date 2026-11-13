"""Funciones para el análisis de expresión diferencial con InMoose."""

from __future__ import annotations

import numpy as np
import pandas as pd

# pandas 3.x entrega cadenas como StringDtype e InMoose no lo interpreta;
# hay que desactivarlo ANTES de importar inmoose
pd.set_option("future.infer_string", False)

from inmoose.deseq2 import DESeqDataSet, nbinomWaldTest

# De acuerdo con Nobori y colaboradores (2018), se utilizaran los siguientes genes
# de referencia: los cuatro primeros deben inducirse in planta y los
# dos últimos mantenerse estables (Nobori et al., 2018)
CONTROLES_POSITIVOS = {
    "PSPTO1404": "hrpL",
    "PSPTO4001": "avrPto",
    "PSPTO4709": "cmaA",
    "PSPTO2133": "pvdS",
}
CONTROLES_NEGATIVOS = {
    "PSPTO1287": "gapA",
    "PSPTO1745": "gyrA",
}


def build_dataset(counts: pd.DataFrame,
                  meta: pd.DataFrame,
                  design: str = "~condicion",
                  reference: str = "culture") -> DESeqDataSet:
    """Construye el objeto del modelo a partir de conteos crudos.

    InMoose espera las muestras en las filas y los genes en las columnas,
    al revés de como se guarda la matriz, por eso se realizará la transposición.

    Parameters
    ----------
    counts : pandas.DataFrame
        Conteos crudos enteros, genes × muestras.
    meta : pandas.DataFrame
        Metadatos indexados por muestra, con la columna ``condition``.
    design : str
        Fórmula del modelo: indica de qué variable depende la expresión,
        en este caso la condición. Usa la notación de patsy porque InMoose 
        internamente usa ``patsy.dmatrix``.
    reference : str
        Nivel de la condición que se usa como denominador del contraste.

    Returns
    -------
    inmoose.deseq2.DESeqDataSet
        Objeto listo para normalizar.

    Raises
    ------
    TypeError
        Si alguna entrada no es DataFrame.
    ValueError
        Si las muestras no coinciden, si hay conteos no enteros o negativos,
        o si falta el nivel de referencia.
    """
    # Validaciones de entrada
    if not isinstance(counts, pd.DataFrame) or not isinstance(meta, pd.DataFrame):
        raise TypeError("counts y meta deben ser DataFrame")
    if set(counts.columns) != set(meta.index):
        raise ValueError(
            "Las muestras de counts y meta no coinciden. Solución:"
            "Asegúrate de que ambas tablas provengan del mismo paso del notebook 01_descarga.ipynb."
        )
    if "condition" not in meta.columns:
        raise ValueError(
            "meta no tiene la columna 'condition'. Solución: "
            "usa el archivo sample_metadata.csv que produce el notebook 01_descarga.ipynb."
        )
    # Verificar que no haya conteos negativos ni decimales, porque InMoose espera conteos enteros
    counts = counts[meta.index]
    if (counts < 0).any().any():
        raise ValueError("Hay conteos negativos en la matriz")
    if not np.allclose(counts.values, np.round(counts.values)):
        raise ValueError(
            "La matriz contiene valores decimales. Solución: usa los conteos "
            "crudos de counts_filtered.csv, no CPM ni valores normalizados."
        )
    # Verifica que únicamente haya un nivel de referencia y que esté presente en la columna condition
    niveles = list(pd.unique(meta["condition"]))
    if reference not in niveles:
        raise ValueError(
            f"El nivel de referencia '{reference}' no está en condition "
            f"{niveles}. Solución: usa uno de esos valores."
        )
    # reordena los niveles, primero el de referencia (denominador de contraste)
    # y luego los demás en el orden que aparecen en la tabla (en este caso solo plant)
    orden = [reference] + [n for n in niveles if n != reference]

    # tabla de muestras, con la condición como categoría y el orden de niveles correcto
    clinical = pd.DataFrame(
        {"condicion": pd.Categorical(meta["condition"].values, categories=orden)},
        index=meta.index,
    )
    # Matriz de conteos, con muestras en filas y genes en columnas, como espera InMoose
    matrix = pd.DataFrame(
        np.array(counts.T.values, dtype=np.int64, copy=True),
        index=meta.index,
        columns=counts.index,
    )
    return DESeqDataSet(countData=matrix, clinicalData=clinical, design=design)


def normalize(dds: DESeqDataSet) -> tuple[DESeqDataSet, pd.Series]:
    """Calcula los factores de tamaño por el método de mediana de razones.

    Parameters
    ----------
    dds : inmoose.deseq2.DESeqDataSet
        Objeto construido con :func:`build_dataset`.

    Returns
    -------
    dds : inmoose.deseq2.DESeqDataSet
        El mismo objeto con los factores de tamaño calculados.
    factors : pandas.Series
        Factor de tamaño de cada muestra.

    Raises
    ------
    RuntimeError
        Si InMoose no pudo calcular los factores de tamaño.
    """
    dds = dds.estimateSizeFactors()
    if "sizeFactor" in dds.obs:
        factors = pd.Series(dds.obs["sizeFactor"], index=dds.obs_names)
    elif "sizeFactors" in dds.obs:
        factors = pd.Series(dds.obs["sizeFactors"], index=dds.obs_names)
    else:
        raise RuntimeError(
            "InMoose no dejó los factores de tamaño en dds.obs. Solución: "
            "verifica la versión instalada con 'pip show inmoose'."
        )
    return dds, factors.rename("factor_de_tamano")


def estimate_dispersions(dds: DESeqDataSet) -> tuple[DESeqDataSet, pd.DataFrame]:
    """Estima la dispersión de cada gen con encogimiento bayesiano.

    Con pocas réplicas la dispersión de un gen no es estimable por sí sola,
    así que se ajusta una tendencia con todos los genes y cada estimación
    individual se acerca a esa curva en proporción a su confiabilidad.

    Parameters
    ----------
    dds : inmoose.deseq2.DESeqDataSet
        Objeto ya normalizado.

    Returns
    -------
    dds : inmoose.deseq2.DESeqDataSet
        El mismo objeto con las dispersiones estimadas.
    table : pandas.DataFrame
        Columnas ``baseMean``, ``dispGeneEst`` (estimación individual),
        ``dispFit`` (valor de la tendencia) y ``dispersion`` (final).
    """
    dds = dds.estimateDispersions()
    columnas = ["baseMean", "dispGeneEst", "dispFit", "dispersion"]
    table = pd.DataFrame(
        {c: np.asarray(dds.var[c]) for c in columnas if c in dds.var},
        index=dds.var_names,
    )
    return dds, table


def fit_model(dds: DESeqDataSet) -> DESeqDataSet:
    """Ajusta el modelo lineal generalizado y aplica la prueba de Wald.

    Parameters
    ----------
    dds : inmoose.deseq2.DESeqDataSet
        Objeto con factores de tamaño y dispersiones ya estimados.

    Returns
    -------
    inmoose.deseq2.DESeqDataSet
        Objeto con los coeficientes ajustados.
    """
    return nbinomWaldTest(dds, quiet=True)


def extract_results(dds: DESeqDataSet,
                    gene_names: pd.Index,
                    contrast: list[str] | None = None,
                    alpha: float = 0.05) -> pd.DataFrame:
    """Extrae la tabla de resultados del contraste de interés.

    Parameters
    ----------
    dds : inmoose.deseq2.DESeqDataSet
        Objeto con el modelo ya ajustado.
    gene_names : pandas.Index
        Identificadores de los genes, en el orden de la matriz.
    contrast : list of str, optional
        ``[variable, numerador, denominador]``. Por omisión compara
        ``plant`` contra ``culture``.
    alpha : float
        Umbral de FDR que usa el filtrado independiente.

    Returns
    -------
    pandas.DataFrame
        Columnas ``baseMean``, ``log2FoldChange``, ``lfcSE``, ``stat``,
        ``pvalue`` y ``adj_pvalue``, indexadas por gen.

    Raises
    ------
    ValueError
        Si alpha está fuera del intervalo (0, 1).
    """
    if not 0 < alpha < 1:
        raise ValueError("alpha debe estar entre 0 y 1")
    if contrast is None:
        contrast = ["condicion", "plant", "culture"]

    results = pd.DataFrame(dds.results(contrast=contrast, alpha=alpha))
    results.index = gene_names
    results.index.name = "locus_key"
    return results


def classify_degs(results: pd.DataFrame,
                  alpha: float = 0.05,
                  lfc_min: float = 1.0) -> pd.DataFrame:
    """Marca cada gen como inducido, reprimido o sin cambio.

    Parameters
    ----------
    results : pandas.DataFrame
        Tabla devuelta por :func:`extract_results`.
    alpha : float
        Umbral de FDR.
    lfc_min : float
        Cambio mínimo en log2 para considerar un gen diferencial.

    Returns
    -------
    pandas.DataFrame
        La misma tabla con la columna ``regulacion``.

    Raises
    ------
    ValueError
        Si faltan las columnas necesarias.
    """
    requeridas = {"log2FoldChange", "adj_pvalue"}
    if not requeridas.issubset(results.columns):
        raise ValueError(
            f"Faltan columnas {requeridas - set(results.columns)}. Solución: "
            "usa la tabla que devuelve extract_results."
        )

    results = results.copy()
    # adj_pvalue es NaN en los genes que el filtrado independiente descartó;
    # comparar NaN siempre da False, así que quedan como 'sin cambio'
    significativo = (results["adj_pvalue"] < alpha) & (
        results["log2FoldChange"].abs() > lfc_min
    )
    results["regulacion"] = np.where(
        significativo & (results["log2FoldChange"] > 0), "inducido",
        np.where(significativo, "reprimido", "sin cambio"),
    )
    return results


def annotate_results(results: pd.DataFrame,
                     annotation: pd.DataFrame,
                     key: str = "locus_key") -> pd.DataFrame:
    """Une los resultados con la anotación funcional del genoma.

    Parameters
    ----------
    results : pandas.DataFrame
        Tabla de resultados indexada por identificador de gen.
    annotation : pandas.DataFrame
        Anotación con la columna indicada en ``key``.
    key : str
        Columna de la anotación que corresponde al índice de ``results``.

    Returns
    -------
    pandas.DataFrame
        Resultados con las columnas ``gene``, ``final_product`` y
        ``replicon``, ordenados por significancia.

    Raises
    ------
    ValueError
        Si la anotación no tiene la columna clave o no cubre ningún gen.
    """
    if key not in annotation.columns:
        raise ValueError(
            f"La anotación no tiene la columna '{key}'. Solución: usa el "
            "final_annotation.csv que produce el notebook 01."
        )

    columnas = [c for c in ("gene", "final_product", "replicon")
                if c in annotation.columns]
    anotacion = annotation.drop_duplicates(key).set_index(key)[columnas]

    unidos = results.join(anotacion, how="left")
    cubiertos = unidos["final_product"].notna().sum() if "final_product" in unidos else 0
    if cubiertos == 0:
        raise ValueError(
            "Ningún gen coincidió con la anotación. Solución: revisa que los "
            "identificadores tengan el mismo formato (PSPTO0001)."
        )
    return unidos.sort_values("adj_pvalue")


def validate_controls(results: pd.DataFrame) -> pd.DataFrame:
    """Comprueba que los genes de referencia se comporten como se espera.

    Es la validación biológica del análisis: si los genes de virulencia no
    aparecen inducidos, el problema está en el diseño o en los datos, no en
    la interpretación.

    Parameters
    ----------
    results : pandas.DataFrame
        Tabla devuelta por :func:`extract_results`.

    Returns
    -------
    pandas.DataFrame
        Una fila por gen control, con su log2FoldChange, su FDR y si cumple
        lo esperado.
    """
    filas = []
    for locus, nombre in {**CONTROLES_POSITIVOS, **CONTROLES_NEGATIVOS}.items():
        esperado = "inducido" if locus in CONTROLES_POSITIVOS else "estable"
        if locus not in results.index:
            filas.append({"gen": nombre, "tipo": esperado, "cumple": "ausente"})
            continue
        lfc = results.loc[locus, "log2FoldChange"]
        padj = results.loc[locus, "adj_pvalue"]
        if esperado == "inducido":
            cumple = bool(lfc > 2 and pd.notna(padj) and padj < 0.05)
        else:
            cumple = bool(abs(lfc) < 2)
        filas.append({
            "gen": nombre,
            "tipo": esperado,
            "log2FoldChange": round(float(lfc), 2),
            "adj_pvalue": float(padj) if pd.notna(padj) else np.nan,
            "cumple": "sí" if cumple else "NO",
        })
    return pd.DataFrame(filas)