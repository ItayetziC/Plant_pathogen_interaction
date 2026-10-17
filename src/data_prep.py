"""Funciones para validar la matriz, seleccionar muestras y limpiar."""

import pandas as pd


def split_id_description(index):
    """Separa cada ID de fila en locus_tag y descripcion (si la trae)."""
    ids = pd.Series(list(index)).astype(str)
    # n=1: corta solo en el primer guion; el resto es la descripcion
    parts = ids.str.split("-", n=1, expand=True)
    table = pd.DataFrame(
        {"original_id": ids, "locus_tag": parts[0].str.strip()}
    )
    # los IDs sin descripcion (p. ej. PSPTOA0002) quedan vacios
    if 1 in parts.columns:
        table["description"] = parts[1].str.strip()
    else:
        table["description"] = None
    # clave sin guion bajo (PSPTO_0001 -> PSPTO0001) para unir con GenBank
    table["locus_key"] = table["locus_tag"].str.replace("_", "", regex=False)
    return table


def validate_gene_ids(ids_table):
    """Comprueba que los locus tag tengan forma PSPTO0001 o PSPTOA0002."""
    keys = ids_table["locus_key"]
    # patron: PSPTO + A o B opcional (plasmidos) + solo digitos
    valid = keys.str.match(r"^PSPTO[AB]?\d+$")
    n_chromosome = int(keys.str.match(r"^PSPTO\d+$").sum())
    n_plasmid_a = int(keys.str.startswith("PSPTOA").sum())
    n_plasmid_b = int(keys.str.startswith("PSPTOB").sum())
    n_described = int(ids_table["description"].notna().sum())
    print(f"Locus tag con formato PSPTO: {int(valid.sum())} de {len(keys)}")
    print(
        f"Cromosoma: {n_chromosome} | Plásmido A: {n_plasmid_a} "
        f"| Plásmido B: {n_plasmid_b}"
    )
    print(f"Genes con descripción: {n_described}")
    print(f"Locus tag repetidos: {int(keys.duplicated().sum())}")
    # si alguno no cumple, muestra ejemplos
    if not valid.all():
        examples = list(ids_table.loc[~valid, "locus_tag"].head(10))
        print("Locus tag con otro formato (ejemplos):", examples)
    return valid


def classify_sample(name):
    """Clasifica una muestra por su nombre: plant, culture o exclude."""
    # Pto_col_A26 -> ["Pto", "col", "A26"]
    parts = name.split("_")
    # solo cepa silvestre (Pto-AvrRpt2 lleva guion y no pasa)
    if parts[0] != "Pto":
        return "exclude"
    # Col-0 sin pretratamiento: exactamente 3 partes
    if len(parts) == 3 and parts[1].lower() == "col":
        return "plant"
    # cultivo en medio King's B (KB en el nombre): exactamente 3 partes
    if len(parts) == 3 and parts[1].lower() == "kb":
        return "culture"
    # mock, flg22, SA, chitin, mutantes, MM, etc.
    return "exclude"


def select_samples(columns):
    """Devuelve una tabla (sample, condition) con cada columna clasificada."""
    table = pd.DataFrame({"sample": list(columns)})
    table["condition"] = table["sample"].apply(classify_sample)
    return table


def verify_selection(selection, expected):
    """Comprueba el conteo por condicion contra el Dataset S1."""
    counts = selection["condition"].value_counts()
    for condition, expected_n in expected.items():
        real_n = int(counts.get(condition, 0))
        # no coincide: detener con mensaje claro
        if real_n != expected_n:
            raise ValueError(
                f"{condition}: se esperaban {expected_n} y hay {real_n}. "
                "Revisa classify_sample() contra los nombres reales."
            )
    print("Selección verificada:", {c: int(counts[c]) for c in expected})


def add_gsm(meta, samples):
    """Añade el codigo GSM a cada muestra (ignorando mayusculas)."""
    # diccionario: titulo en minusculas -> codigo GSM
    lookup = samples.assign(key=samples["title"].str.lower())
    lookup = lookup.set_index("key")["gsm"]
    meta = meta.copy()
    meta["gsm"] = meta.index.str.lower().map(lookup)
    # si alguna no se encontro, detener
    if meta["gsm"].isna().any():
        missing = list(meta.index[meta["gsm"].isna()])
        raise ValueError(f"Hay muestras sin GSM: {missing}")
    return meta


def clean_and_filter(
    matrix, ids_table, selection, min_cpm=1, min_samples=2
):
    """Deja las muestras elegidas, limpia IDs y quita genes sin expresion."""
    chosen = selection[selection["condition"] != "exclude"]
    counts = matrix[chosen["sample"].tolist()].copy()
    # reemplaza el ID largo por el locus tag limpio (mismo orden de filas)
    counts.index = ids_table["locus_key"].values
    if counts.index.duplicated().any():
        raise ValueError("Hay locus tag duplicados tras la limpieza")
    if counts.isna().any().any():
        raise ValueError("Hay NaN en los conteos")
    print(
        f"Muestras: {counts.shape[1]} | "
        f"Genes antes del filtro: {counts.shape[0]}"
    )

    library_size = counts.sum()
    print("\nTamaño de librería (conteos totales por muestra):")
    print(library_size.sort_values().to_string())

    # CPM: conteos por millon, corrige la profundidad de cada muestra
    cpm = counts / library_size * 1e6
    # un gen se queda con CPM > min_cpm en al menos min_samples muestras
    keep = (cpm > min_cpm).sum(axis=1) >= min_samples
    filtered = counts[keep]
    print(f"\nGenes que pasan el filtro: {int(keep.sum())} de {len(keep)}")

    # condicion de cada muestra elegida, en el mismo orden que las columnas
    meta = chosen.set_index("sample")[["condition"]]
    return filtered, meta