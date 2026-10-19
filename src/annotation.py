"""Funciones para descargar la anotación de DC3000 y unirla a la matriz."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import pandas as pd
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

from geo_download import entrez_query


def _read_genbank(path: str | Path) -> SeqRecord:
    """Lee un archivo GenBank con un solo registro.

    Parameters
    ----------
    path : str or pathlib.Path
        Ruta del archivo ``.gb``.

    Returns
    -------
    Bio.SeqRecord.SeqRecord
        Registro leído.

    Raises
    ------
    FileNotFoundError
        Si el archivo no existe.
    ValueError
        Si el archivo no tiene formato GenBank válido.
    """
    try:
        return SeqIO.read(path, "genbank")
    except FileNotFoundError as error:
        raise FileNotFoundError(
            f"No existe: {path}. Solución: ejecuta antes download_genbank()."
        ) from error
    except ValueError as error:
        raise ValueError(
            f"{path} no es un GenBank válido ({error}). Solución: bórralo "
            "y vuelve a descargarlo."
        ) from error


def validate_nuccore_accessions(accessions: list[str]) -> list[str]:
    """Comprueba que cada acceso exista en NCBI y sea de DC3000.

    Parameters
    ----------
    accessions : list of str
        Accesos de nuccore, con o sin versión.

    Returns
    -------
    list of str
        Accesos con la versión vigente, p. ej. ``"NC_004578.1"``.

    Raises
    ------
    TypeError
        Si ``accessions`` no es una lista de texto.
    ValueError
        Si falta algún registro o alguno no es de DC3000.
    """
    if not isinstance(accessions, list) or not all(
        isinstance(acc, str) for acc in accessions
    ):
        raise TypeError("accessions debe ser una lista de str")
    # una sola consulta para todos los accesos reduce llamadas a NCBI
    summaries = entrez_query(
        Entrez.esummary, db="nuccore", id=",".join(accessions)
    )
    if len(summaries) != len(accessions):
        raise ValueError(
            f"Se pidieron {len(accessions)} accesos y NCBI devolvió "
            f"{len(summaries)}"
        )
    valid = []
    for record in summaries:
        title = str(record["Title"])
        version = str(record["AccessionVersion"])
        print(version, "|", title[:80])
        # un acceso mal escrito puede existir pero ser de otro organismo
        if "DC3000" not in title:
            raise ValueError(f"{version} no parece ser DC3000: {title}")
        valid.append(version)
    return valid


def download_genbank(accession: str, folder: str | Path) -> Path:
    """Descarga un registro GenBank y comprueba que tenga CDS.

    Parameters
    ----------
    accession : str
        Acceso de nuccore con versión.
    folder : str or pathlib.Path
        Carpeta destino; se crea si no existe.

    Returns
    -------
    pathlib.Path
        Ruta del archivo ``.gb``.

    Raises
    ------
    ValueError
        Si el archivo no es GenBank válido o no contiene CDS.
    """
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    destination = folder / f"{accession}.gb"
    # evita repetir una consulta pesada a NCBI si ya se tiene el archivo
    if destination.exists() and destination.stat().st_size > 0:
        print(f"Ya existe: {destination}")
    else:
        # gbwithparts incluye las CDS; "gb" a secas puede omitirlas
        text = entrez_query(
            Entrez.efetch, parse=False, db="nuccore", id=accession,
            rettype="gbwithparts", retmode="text",
        )
        destination.write_text(text, encoding="utf-8")
        print(f"Descargado: {destination}")
    record = _read_genbank(destination)
    n_cds = sum(1 for feature in record.features if feature.type == "CDS")
    if n_cds == 0:
        raise ValueError(
            f"{destination} no contiene CDS. Solución: bórralo y vuelve a "
            "descargarlo."
        )
    print(f"  {record.id}: {len(record.seq)} pb, {n_cds} CDS")
    return destination


def extract_cds(paths: Iterable[str | Path]) -> pd.DataFrame:
    """Lee archivos GenBank y devuelve una fila por CDS.

    Parameters
    ----------
    paths : iterable of str or pathlib.Path
        Rutas de los archivos ``.gb``.

    Returns
    -------
    pandas.DataFrame
        Columnas ``locus_tag``, ``old_locus_tag``, ``gene``, ``product``,
        ``note``, ``replicon``, ``start``, ``end`` y ``strand``.

    Raises
    ------
    ValueError
        Si no se encontró ninguna CDS.
    """
    rows = []
    for path in paths:
        record = _read_genbank(path)
        for feature in record.features:
            if feature.type != "CDS":
                continue
            quals = feature.qualifiers
            rows.append(
                {
                    "locus_tag": quals.get("locus_tag", [None])[0],
                    "old_locus_tag": quals.get("old_locus_tag", [None])[0],
                    "gene": quals.get("gene", [None])[0],
                    "product": quals.get("product", [None])[0],
                    "note": "; ".join(quals.get("note", [])),
                    "replicon": record.id,
                    # Biopython cuenta desde 0; GenBank y NCBI desde 1
                    "start": int(feature.location.start) + 1,
                    "end": int(feature.location.end),
                    "strand": feature.location.strand,
                }
            )
    if not rows:
        raise ValueError("No se encontraron CDS en los archivos dados")
    return pd.DataFrame(rows)


def merge_with_matrix(annotation: pd.DataFrame,
                      ids_table: pd.DataFrame) -> pd.DataFrame:
    """Une GenBank con la matriz por ``old_locus_tag`` o ``locus_tag``.

    Parameters
    ----------
    annotation : pandas.DataFrame
        Salida de :func:`extract_cds`.
    ids_table : pandas.DataFrame
        Tabla con la columna ``locus_key``.

    Returns
    -------
    pandas.DataFrame
        Anotación de los genes con pareja, indexada por ``locus_key``.

    Raises
    ------
    TypeError
        Si alguna entrada no es DataFrame.
    ValueError
        Si faltan columnas necesarias.
    """
    for table, name in ((annotation, "annotation"), (ids_table, "ids_table")):
        if not isinstance(table, pd.DataFrame):
            raise TypeError(f"{name} debe ser DataFrame")
    if "locus_key" not in ids_table.columns:
        raise ValueError("A ids_table le falta la columna locus_key")
    annotation = annotation.copy()
    for column in ("locus_tag", "old_locus_tag"):
        annotation["key_" + column] = annotation[column].str.replace(
            "_", "", regex=False
        )
    # la matriz usa el codigo antiguo; locus_tag solo es respaldo
    mapping = {}
    for column in ("key_old_locus_tag", "key_locus_tag"):
        for key, row in zip(annotation[column], annotation.index):
            if isinstance(key, str):
                mapping.setdefault(key, row)
    positions = ids_table["locus_key"].map(mapping)
    n_found = int(positions.notna().sum())
    print(
        f"Genes de la matriz con anotación en GenBank: {n_found} "
        f"de {len(ids_table)}"
    )
    missing = ids_table.loc[positions.isna(), "locus_key"]
    if len(missing) > 0:
        print("Sin anotación (ejemplos):", list(missing.head(10)))
    columns = [
        "locus_tag", "old_locus_tag", "gene", "product", "note",
        "replicon", "start", "end", "strand",
    ]
    merged = annotation.loc[positions.dropna().astype(int), columns].copy()
    merged.index = ids_table.loc[positions.notna(), "locus_key"].values
    merged.index.name = "locus_key"
    return merged


def combine_annotation(annotation: pd.DataFrame,
                       ids_table: pd.DataFrame) -> pd.DataFrame:
    """Combina GenBank con la descripción de la matriz, una fila por gen.

    Parameters
    ----------
    annotation : pandas.DataFrame
        Salida de :func:`merge_with_matrix`.
    ids_table : pandas.DataFrame
        Tabla con ``locus_key``, ``locus_tag`` y ``description``.

    Returns
    -------
    pandas.DataFrame
        Una fila por gen con ``final_product`` (GenBank si existe; si
        no, la matriz) y ``source``.

    Raises
    ------
    ValueError
        Si faltan columnas o la unión cambia el número de genes.
    """
    needed = ["locus_key", "locus_tag", "description"]
    if not set(needed) <= set(ids_table.columns):
        raise ValueError(f"ids_table debe tener las columnas {needed}")
    base = ids_table.set_index("locus_key")[["locus_tag", "description"]]
    # sufijos de origen para no confundir las dos fuentes en el informe
    base = base.rename(
        columns={
            "locus_tag": "locus_tag_matrix",
            "description": "description_matrix",
        }
    )
    genbank = annotation.rename(
        columns={
            "locus_tag": "locus_tag_refseq",
            "product": "product_genbank",
            "note": "note_genbank",
        }
    )
    columns = [
        "locus_tag_refseq", "gene", "product_genbank", "note_genbank",
        "replicon", "start", "end", "strand",
    ]
    # left: los genes sin pareja en GenBank se conservan con NaN
    table = base.join(genbank[columns], how="left")
    if len(table) != len(ids_table):
        raise ValueError("La unión cambió el número de genes")
    # GenBank primero porque RefSeq es la anotación más reciente
    table["final_product"] = table["product_genbank"].fillna(
        table["description_matrix"]
    )
    has_genbank = table["product_genbank"].notna()
    table["source"] = has_genbank.map(
        {True: "GenBank (RefSeq)", False: "Count matrix"}
    )
    both = has_genbank & table["description_matrix"].notna()
    same = (
        table.loc[both, "product_genbank"].str.lower().str.strip()
        == table.loc[both, "description_matrix"].str.lower().str.strip()
    )
    print(
        f"Genes con ambas fuentes: {int(both.sum())} "
        f"| descripción idéntica: {int(same.sum())}"
    )
    sources = table["source"].value_counts().to_dict()
    print(f"Fuente del producto final: {sources}")
    n_empty = int(table["final_product"].isna().sum())
    print(f"Genes sin descripción en ninguna fuente: {n_empty}")
    return table
