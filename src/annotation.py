"""Funciones para descargar la anotacion de DC3000 y unirla a la matriz."""

from contextlib import closing
from pathlib import Path

import pandas as pd
from Bio import Entrez, SeqIO


def validate_nuccore_accessions(accessions):
    """Comprueba que cada acceso exista en NCBI y sea de DC3000."""
    # esummary pide la ficha de varios registros a la vez
    with closing(
        Entrez.esummary(db="nuccore", id=",".join(accessions))
    ) as handle:
        summaries = Entrez.read(handle)
    # NCBI debe devolver una ficha por acceso
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
        # debe ser de la cepa correcta
        if "DC3000" not in title:
            raise ValueError(f"{version} no parece ser DC3000: {title}")
        valid.append(version)
    return valid


def download_genbank(accession, folder):
    """Descarga un registro GenBank y comprueba que el formato sea valido."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    destination = folder / f"{accession}.gb"
    # si ya existe, no se descarga otra vez
    if destination.exists() and destination.stat().st_size > 0:
        print(f"Ya existe: {destination}")
    else:
        # efetch descarga el registro completo en formato GenBank
        with closing(
            Entrez.efetch(
                db="nuccore",
                id=accession,
                rettype="gbwithparts",
                retmode="text",
            )
        ) as handle:
            text = handle.read()
        destination.write_text(text, encoding="utf-8")
        print(f"Descargado: {destination}")
    # si el formato no es GenBank valido, SeqIO lanza error
    record = SeqIO.read(destination, "genbank")
    n_cds = sum(1 for feature in record.features if feature.type == "CDS")
    # un genoma sin CDS indica una descarga incorrecta
    if n_cds == 0:
        raise ValueError(f"{destination} no contiene CDS")
    print(f"  {record.id}: {len(record.seq)} pb, {n_cds} CDS")
    return destination


def extract_cds(paths):
    """Lee los archivos GenBank y devuelve una fila por CDS."""
    rows = []
    for path in paths:
        record = SeqIO.read(path, "genbank")
        for feature in record.features:
            # solo interesan los genes codificantes
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
                    # inicio desde 1, como en los archivos GenBank
                    "start": int(feature.location.start) + 1,
                    "end": int(feature.location.end),
                    # 1 = hebra directa, -1 = reversa
                    "strand": feature.location.strand,
                }
            )
    return pd.DataFrame(rows)


def merge_with_matrix(annotation, ids_table):
    """Une GenBank con la matriz (prueba old_locus_tag y locus_tag)."""
    annotation = annotation.copy()
    # claves sin guion bajo, igual que en la matriz
    for column in ("locus_tag", "old_locus_tag"):
        annotation["key_" + column] = annotation[column].str.replace(
            "_", "", regex=False
        )
    # diccionario: clave -> fila de la anotacion (primero el codigo antiguo)
    mapping = {}
    for column in ("key_old_locus_tag", "key_locus_tag"):
        for key, row in zip(annotation[column], annotation.index):
            if isinstance(key, str):
                mapping.setdefault(key, row)
    # fila de anotacion de cada gen (NaN si no hay)
    positions = ids_table["locus_key"].map(mapping)
    n_found = int(positions.notna().sum())
    print(
        f"Genes de la matriz con anotación en GenBank: {n_found} "
        f"de {len(ids_table)}"
    )
    missing = ids_table.loc[positions.isna(), "locus_key"]
    # si hay, muestra ejemplos para diagnosticar
    if len(missing) > 0:
        print("Sin anotación (ejemplos):", list(missing.head(10)))
    columns = [
        "locus_tag", "old_locus_tag", "gene", "product", "note",
        "replicon", "start", "end", "strand",
    ]
    merged = annotation.loc[positions.dropna().astype(int), columns].copy()
    # indexa por locus_key de la matriz
    merged.index = ids_table.loc[positions.notna(), "locus_key"].values
    merged.index.name = "locus_key"
    return merged


def combine_annotation(annotation, ids_table):
    """Une GenBank con la descripcion de la matriz, una fila por gen."""
    # base: un gen por fila, con la descripcion que traia la matriz
    base = ids_table.set_index("locus_key")[["locus_tag", "description"]]
    base = base.rename(
        columns={
            "locus_tag": "locus_tag_matrix",
            "description": "description_matrix",
        }
    )
    # nombres de columna que dejan claro de donde viene cada dato
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
    # sin pareja en GenBank queda vacio (NaN)
    table = base.join(genbank[columns], how="left")
    # la union no debe crear ni perder genes
    if len(table) != len(ids_table):
        raise ValueError("La unión cambió el número de genes")
    # producto final: GenBank si existe; si no, la matriz
    table["final_product"] = table["product_genbank"].fillna(
        table["description_matrix"]
    )
    # fuente del producto final, para el informe
    has_genbank = table["product_genbank"].notna()
    table["source"] = has_genbank.map(
        {True: "GenBank (RefSeq)", False: "Count matrix"}
    )
    # concordancia entre fuentes en los genes que tienen ambas
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