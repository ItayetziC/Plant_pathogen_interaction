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
