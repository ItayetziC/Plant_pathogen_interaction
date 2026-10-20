""""Funciones para validar la matriz, elegir muestras, limpiar y filtrar."""

from __future__ import annotations

import warnings
from collections.abc import Iterable

import pandas as pd


def _require_columns(table: pd.DataFrame, columns: Iterable[str],
                     name: str) -> None:
    """Comprueba que una tabla sea DataFrame y tenga ciertas columnas.

    Parameters
    ----------
    table : pandas.DataFrame
        Tabla a revisar.
    columns : iterable of str
        Columnas obligatorias.
    name : str
        Nombre del argumento, para el mensaje de error.

    Raises
    ------
    TypeError
        Si ``table`` no es un DataFrame.
    ValueError
        Si falta alguna columna.
    """
    if not isinstance(table, pd.DataFrame):
        raise TypeError(
            f"{name} debe ser DataFrame, no {type(table).__name__}"
        )
    missing = [column for column in columns if column not in table.columns]
    if missing:
        raise ValueError(
            f"A {name} le faltan las columnas {missing}. Solución: usa la "
            "tabla que genera el paso anterior del notebook."
        )


def split_id_description(index: Iterable[str]) -> pd.DataFrame:
    """Separa cada ID de fila en locus tag y descripción.

    Parameters
    ----------
    index : iterable of str
        IDs de la matriz, p. ej. ``"PSPTO_0001-chromosomal ..."``.

    Returns
    -------
    pandas.DataFrame
        Columnas ``original_id``, ``locus_tag``, ``description`` (None
        si no trae) y ``locus_key`` (locus tag sin guion bajo).
    """
    ids = pd.Series(list(index)).astype(str)
    # n=1 porque las descripciones tambien pueden contener guiones
    parts = ids.str.split("-", n=1, expand=True)
    table = pd.DataFrame(
        {"original_id": ids, "locus_tag": parts[0].str.strip()}
    )
    # si ningun ID trae descripcion, split no crea la segunda columna
    if 1 in parts.columns:
        table["description"] = parts[1].str.strip()
    else:
        table["description"] = None
    # GenBank escribe el old_locus_tag sin guion bajo; asi coinciden
    table["locus_key"] = table["locus_tag"].str.replace("_", "", regex=False)
    return table


def validate_gene_ids(ids_table: pd.DataFrame) -> pd.Series:
    """Comprueba que los locus tag tengan forma PSPTO0001 o PSPTOA0002.

    Parameters
    ----------
    ids_table : pandas.DataFrame
        Salida de :func:`split_id_description`.

    Returns
    -------
    pandas.Series of bool
        True para cada locus tag con formato válido.

    Warns
    -----
    UserWarning
        Si hay locus tag con otro formato o repetidos.
    """
    _require_columns(ids_table, ["locus_key", "locus_tag", "description"],
                     "ids_table")
    keys = ids_table["locus_key"]
    # A y B son los plasmidos; el cromosoma no lleva letra
    valid = keys.str.match(r"^PSPTO[AB]?\d+$")
    n_chromosome = int(keys.str.match(r"^PSPTO\d+$").sum())
    n_plasmid_a = int(keys.str.startswith("PSPTOA").sum())
    n_plasmid_b = int(keys.str.startswith("PSPTOB").sum())
    n_described = int(ids_table["description"].notna().sum())
    n_repeated = int(keys.duplicated().sum())
    print(f"Locus tag con formato PSPTO: {int(valid.sum())} de {len(keys)}")
    print(
        f"Cromosoma: {n_chromosome} | Plásmido A: {n_plasmid_a} "
        f"| Plásmido B: {n_plasmid_b}"
    )
    print(f"Genes con descripción: {n_described}")
    print(f"Locus tag repetidos: {n_repeated}")
    if not valid.all():
        examples = list(ids_table.loc[~valid, "locus_tag"].head(10))
        warnings.warn(f"Locus tag con otro formato: {examples}", UserWarning)
    if n_repeated:
        warnings.warn("Hay locus tag repetidos; la unión con GenBank "
                      "fallará.", UserWarning)
    return valid


def classify_sample(name: str) -> str:
    """Clasifica una muestra por su nombre.

    Parameters
    ----------
    name : str
        Nombre de la columna, p. ej. ``"Pto_col_A26"``.

    Returns
    -------
    str
        ``"plant"`` (Col-0 sin pretratamiento), ``"culture"`` (medio
        King's B) o ``"exclude"``.

    Raises
    ------
    TypeError
        Si ``name`` no es texto.
    """
    if not isinstance(name, str):
        raise TypeError(f"name debe ser str, no {type(name).__name__}")
    parts = name.split("_")
    # Pto-AvrRpt2 y otras cepas usan guion, asi que no pasan este filtro
    if parts[0] != "Pto":
        return "exclude"
    # 3 partes deja fuera pretratamientos (flg22, SA...) y mutantes
    if len(parts) == 3 and parts[1].lower() == "col":
        return "plant"
    if len(parts) == 3 and parts[1].lower() == "kb":
        return "culture"
    return "exclude"


def select_samples(columns: Iterable[str]) -> pd.DataFrame:
    """Clasifica cada columna de la matriz.

    Parameters
    ----------
    columns : iterable of str
        Nombres de las muestras.

    Returns
    -------
    pandas.DataFrame
        Columnas ``sample`` y ``condition``.
    """
    table = pd.DataFrame({"sample": list(columns)})
    table["condition"] = table["sample"].apply(classify_sample)
    return table


def verify_selection(selection: pd.DataFrame,
                     expected: dict[str, int]) -> None:
    """Compara el número de muestras por condición con el Dataset S1.

    Parameters
    ----------
    selection : pandas.DataFrame
        Salida de :func:`select_samples`.
    expected : dict of {str: int}
        Muestras esperadas por condición, p. ej. ``{"plant": 13}``.

    Raises
    ------
    ValueError
        Si alguna condición no tiene el número esperado.
    """
    _require_columns(selection, ["condition"], "selection")
    counts = selection["condition"].value_counts()
    for condition, expected_n in expected.items():
        real_n = int(counts.get(condition, 0))
        if real_n != expected_n:
            raise ValueError(
                f"{condition}: se esperaban {expected_n} y hay {real_n}. "
                "Solución: revisa classify_sample() contra los nombres "
                "reales."
            )
    print("Selección verificada:", {c: int(counts[c]) for c in expected})


def add_gsm(meta: pd.DataFrame, samples: pd.DataFrame) -> pd.DataFrame:
    """Añade el código GSM a cada muestra.

    Parameters
    ----------
    meta : pandas.DataFrame
        Metadatos indexados por nombre de muestra.
    samples : pandas.DataFrame
        Tabla de GEO con columnas ``gsm`` y ``title``.

    Returns
    -------
    pandas.DataFrame
        Copia de ``meta`` con la columna ``gsm``.

    Raises
    ------
    ValueError
        Si alguna muestra no tiene GSM.
    """
    _require_columns(samples, ["gsm", "title"], "samples")
    # GEO y la matriz no escriben igual las mayusculas de los titulos
    lookup = samples.assign(key=samples["title"].str.lower())
    lookup = lookup.set_index("key")["gsm"]
    meta = meta.copy()
    meta["gsm"] = meta.index.str.lower().map(lookup)
    if meta["gsm"].isna().any():
        missing = list(meta.index[meta["gsm"].isna()])
        raise ValueError(
            f"Hay muestras sin GSM: {missing}. Solución: compara sus "
            "nombres con la columna title de GEO_samples.csv."
        )
    return meta


def build_count_table(
    matrix: pd.DataFrame,
    ids_table: pd.DataFrame,
    selection: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Deja las muestras elegidas y pone ``locus_key`` como índice.

    No quita genes: la matriz de salida es la que se revisa con
    :func:`quality_control.library_metrics` antes de filtrar.

    Parameters
    ----------
    matrix : pandas.DataFrame
        Matriz cruda genes × muestras.
    ids_table : pandas.DataFrame
        Salida de :func:`split_id_description`, en el orden de ``matrix``.
    selection : pandas.DataFrame
        Salida de :func:`select_samples`.

    Returns
    -------
    counts : pandas.DataFrame
        Conteos sin filtrar de las muestras elegidas, indexados por
        ``locus_key``.
    meta : pandas.DataFrame
        Columna ``condition``, en el orden de las columnas de ``counts``.

    Raises
    ------
    TypeError
        Si ``matrix`` no es DataFrame.
    ValueError
        Si las tablas no concuerdan o hay NaN o duplicados.
    """
    _require_columns(ids_table, ["locus_key"], "ids_table")
    _require_columns(selection, ["sample", "condition"], "selection")
    if not isinstance(matrix, pd.DataFrame):
        raise TypeError("matrix debe ser DataFrame")
    # el indice se reemplaza por posicion, asi que el tamaño debe coincidir
    if len(ids_table) != len(matrix):
        raise ValueError("ids_table y matrix no tienen el mismo número de "
                         "genes. Solución: genera ids_table desde matrix.")
    chosen = selection[selection["condition"] != "exclude"]
    try:
        counts = matrix[chosen["sample"].tolist()].copy()
    except KeyError as error:
        raise ValueError(f"Muestras que no están en la matriz: {error}"
                         ) from error
    counts.index = ids_table["locus_key"].values
    if counts.index.duplicated().any():
        raise ValueError("Hay locus tag duplicados tras la limpieza")
    if counts.isna().any().any():
        raise ValueError("Hay NaN en los conteos")
    print(
        f"Muestras: {counts.shape[1]} | "
        f"Genes antes del filtro: {counts.shape[0]}"
    )
    meta = chosen.set_index("sample")[["condition"]]
    return counts, meta


def filter_low_counts(counts: pd.DataFrame, min_cpm: float = 1,
                      min_samples: int = 2) -> pd.DataFrame:
    """Quita los genes sin expresión suficiente.

    Parameters
    ----------
    counts : pandas.DataFrame
        Conteos genes × muestras, salida de :func:`build_count_table`.
    min_cpm : float, optional
        CPM mínimo para considerar un gen expresado. Por defecto 1.
    min_samples : int, optional
        Muestras que deben superar ``min_cpm``. Por defecto 2, el tamaño
        del grupo más chico.

    Returns
    -------
    pandas.DataFrame
        Conteos de los genes que pasan el filtro.

    Raises
    ------
    TypeError
        Si ``counts`` no es DataFrame o los umbrales no son números.
    ValueError
        Si ``min_samples`` no cabe en el número de muestras o alguna
        muestra no tiene conteos.
    """
    if not isinstance(counts, pd.DataFrame):
        raise TypeError("counts debe ser DataFrame")
    valid_cpm = isinstance(min_cpm, (int, float))
    if not valid_cpm or not isinstance(min_samples, int):
        raise TypeError("min_cpm debe ser número y min_samples entero")
    # con mas muestras exigidas que columnas, ningun gen pasaria
    if not 1 <= min_samples <= counts.shape[1]:
        raise ValueError(
            f"min_samples debe estar entre 1 y {counts.shape[1]}. "
            "Solución: usa el tamaño del grupo más chico (2)."
        )
    library_size = counts.sum()
    # una muestra sin conteos daria division entre cero en el CPM
    if (library_size == 0).any():
        empty = list(library_size.index[library_size == 0])
        raise ValueError(f"Muestras con 0 conteos totales: {empty}")
    # CPM y no conteos: las librerias varian ~4 veces de tamaño
    cpm = counts / library_size * 1e6
    keep = (cpm > min_cpm).sum(axis=1) >= min_samples
    print(f"Genes que pasan el filtro: {int(keep.sum())} de {len(keep)}")
    return counts[keep]
