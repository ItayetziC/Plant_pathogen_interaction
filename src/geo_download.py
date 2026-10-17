"""Funciones para consultar GEO con Entrez y descargar archivos."""

import re
import urllib.error
import urllib.request
from contextlib import closing
from pathlib import Path

import pandas as pd
from Bio import Entrez


def get_gse_summary(accession):
    """Busca una serie GSE en GEO y devuelve su ficha resumen."""
    term = f"{accession}[ACCN] AND gse[ETYP]"
    # esearch devuelve un handle (conexion abierta); closing lo cierra
    with closing(Entrez.esearch(db="gds", term=term)) as handle:
        search = Entrez.read(handle)
    ids = search["IdList"]
    # lista vacia = la serie no existe
    if not ids:
        raise ValueError(f"No se encontró {accession} en GEO")
    # segunda consulta: la ficha resumen del primer ID
    with closing(Entrez.esummary(db="gds", id=ids[0])) as handle:
        summary = Entrez.read(handle)[0]
    # comprobar que sea justo la serie pedida
    found = summary.get("Accession")
    if found != accession:
        raise ValueError(f"Se pidió {accession} pero NCBI devolvió {found}")
    return summary


def validate_bioproject(summary):
    """Busca el BioProject enlazado a la serie y devuelve sus codigos PRJ."""
    # elink: que registros de bioproject estan enlazados a esta serie
    with closing(
        Entrez.elink(dbfrom="gds", db="bioproject", id=summary["Id"])
    ) as handle:
        links = Entrez.read(handle)
    ids = [
        link["Id"]
        for link_set in links
        for database in link_set.get("LinkSetDb", [])
        for link in database["Link"]
    ]
    # sin enlaces: avisa y continua
    if not ids:
        print("Aviso: GEO no enlaza un BioProject a esta serie vía Entrez")
        return {}
    # esummary convierte los IDs internos en codigos PRJ legibles
    with closing(
        Entrez.esummary(db="bioproject", id=",".join(ids))
    ) as handle:
        result = Entrez.read(handle)
    try:
        docs = result["DocumentSummarySet"]["DocumentSummary"]
        return {doc["Project_Acc"]: True for doc in docs}
    except (KeyError, TypeError):
        # si la estructura es otra, no rompe el pipeline
        print("Aviso: no pude leer los códigos PRJ; recibí:", result)
        return {}


def get_samples_table(summary):
    """Convierte la lista de muestras de la ficha en una tabla."""
    rows = [
        {"gsm": sample["Accession"], "title": sample["Title"]}
        for sample in summary["Samples"]
    ]
    # la ficha reporta cuantas muestras debe haber
    expected = int(summary["n_samples"])
    if len(rows) != expected:
        raise ValueError(f"Se esperaban {expected} muestras y hay {len(rows)}")
    return pd.DataFrame(rows)


def list_supplementary_files(summary):
    """Devuelve la URL de la carpeta suppl/ y los archivos que contiene."""
    # misma ruta del FTP, pero abrible como pagina web
    ftp = summary["FTPLink"].replace("ftp://", "https://")
    suppl_url = ftp.rstrip("/") + "/suppl/"
    try:
        with urllib.request.urlopen(suppl_url, timeout=60) as response:
            html = response.read().decode()
    except urllib.error.HTTPError as error:
        raise FileNotFoundError(
            f"No se pudo abrir {suppl_url} (HTTP {error.code})"
        )
    links = re.findall(r'href="([^"?/][^"]*)"', html)
    # un archivo real no contiene ':' ni '/' (descarta enlaces externos)
    files = [link for link in links if ":" not in link and "/" not in link]
    return suppl_url, files


def remote_exists(url):
    """Comprueba si un archivo existe en el servidor sin descargarlo."""
    # HEAD = solo pregunta, no descarga
    request = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(request, timeout=60):
            return True
    except urllib.error.HTTPError:
        return False


def is_gzip(path):
    """Comprueba que un archivo sea gzip mirando sus dos primeros bytes."""
    with open(path, "rb") as file:
        return file.read(2) == b"\x1f\x8b"


def download_file(suppl_url, name, dest_folder):
    """Descarga un archivo si existe en el servidor y aun no lo tienes."""
    folder = Path(dest_folder)
    folder.mkdir(parents=True, exist_ok=True)
    destination = folder / name
    # si ya lo tienes, no se descarga otra vez
    if destination.exists() and destination.stat().st_size > 0:
        print(f"Ya existe: {destination}")
        return destination
    url = suppl_url + name
    # validar antes de descargar
    if not remote_exists(url):
        raise FileNotFoundError(f"No existe en el servidor: {url}")
    urllib.request.urlretrieve(url, destination)
    # validar que no quedo vacio
    size = destination.stat().st_size
    if size == 0:
        raise IOError(f"La descarga quedó vacía: {destination}")
    # validar el formato del archivo
    if name.endswith(".gz") and not is_gzip(destination):
        raise IOError(f"{destination} dice ser .gz pero no es un gzip válido")
    print(f"Descargado: {destination} ({size / 1e6:.2f} MB)")
    return destination


def validate_count_matrix(path):
    """Lee la matriz y comprueba que tenga formato de conteos crudos."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"No existe: {path}")
    # sep tab: pandas descomprime .gz automaticamente
    matrix = pd.read_csv(path, sep="\t", index_col=0)
    # con una sola columna, el separador no era tab
    if matrix.shape[1] < 2:
        raise ValueError("La matriz tiene <2 columnas: revisa el separador")
    numeric = matrix.select_dtypes("number")
    if numeric.shape[1] != matrix.shape[1]:
        raise ValueError("Hay columnas no numéricas en la matriz")
    # sin decimales = conteos crudos
    all_integers = bool(((numeric.fillna(0) % 1) == 0).all().all())
    print(f"Genes x muestras: {matrix.shape}")
    print(f"¿Todos enteros? {all_integers}")
    if not all_integers:
        print("ADVERTENCIA: no parecen conteos crudos; revisa antes de seguir")
    n_nan = int(matrix.isna().sum().sum())
    n_duplicated = int(matrix.index.duplicated().sum())
    print(f"NaN: {n_nan} | IDs duplicados: {n_duplicated}")
    return matrix