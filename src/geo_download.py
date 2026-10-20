"""Funciones para consultar GEO con Entrez y descargar archivos."""

from __future__ import annotations

import http.client
import re
import urllib.error
import urllib.request
import warnings
from collections.abc import Callable
from contextlib import closing
from pathlib import Path
from typing import Any

import pandas as pd
from Bio import Entrez


def entrez_query(
    function: Callable[..., Any], parse: bool = True, **params: Any
) -> Any:
    """Ejecuta una consulta de Entrez y traduce sus errores.

    Parameters
    ----------
    function : callable
        Función de ``Bio.Entrez`` a llamar (``esearch``, ``efetch``...).
    parse : bool, optional
        Si es True, interpreta la respuesta con ``Entrez.read``; si es
        False, devuelve el texto crudo. Por defecto True.
    **params
        Argumentos que se pasan a ``function`` (``db``, ``id``...).

    Returns
    -------
    Any
        Respuesta interpretada por ``Entrez.read`` o texto crudo.

    Raises
    ------
    ValueError
        Si no hay correo configurado o NCBI rechaza la consulta.
    ConnectionError
        Si no hay conexión con NCBI.
    """
    # NCBI bloquea consultas anonimas, mejor avisar antes de intentarlo
    if not Entrez.email:
        raise ValueError(
            "Falta Entrez.email. Solución: asigna tu correo antes de "
            "consultar, p. ej. Entrez.email = 'tu@correo.com'."
        )
    try:
        with closing(function(**params)) as handle:
            return Entrez.read(handle) if parse else handle.read()
    except urllib.error.HTTPError as error:
        raise ValueError(
            f"NCBI rechazó la consulta (HTTP {error.code}). Solución: "
            "revisa el acceso o espera unos segundos si fue un 429."
        ) from error
    except (urllib.error.URLError, http.client.HTTPException) as error:
        raise ConnectionError(
            "No hay conexión con NCBI. Solución: revisa tu internet y "
            "vuelve a ejecutar la celda."
        ) from error
    except RuntimeError as error:
        # Entrez.read lanza RuntimeError cuando NCBI responde con un error
        raise ValueError(f"NCBI devolvió un error: {error}") from error


def get_gse_summary(accession: str) -> dict[str, Any]:
    """Busca una serie GSE en GEO y devuelve su ficha resumen.

    Parameters
    ----------
    accession : str
        Código de la serie, p. ej. ``"GSE103441"``.

    Returns
    -------
    dict
        Ficha de ``esummary`` (``Id``, ``title``, ``Samples``...).

    Raises
    ------
    TypeError
        Si ``accession`` no es texto.
    ValueError
        Si la serie no existe o NCBI devuelve otra distinta.
    """
    if not isinstance(accession, str):
        raise TypeError(
            f"accession debe ser str, no {type(accession).__name__}"
        )
    term = f"{accession}[ACCN] AND gse[ETYP]"
    search = entrez_query(Entrez.esearch, db="gds", term=term)
    ids = search["IdList"]
    if not ids:
        raise ValueError(
            f"No se encontró {accession} en GEO. Solución: revisa que el "
            "código esté bien escrito."
        )
    summary = entrez_query(Entrez.esummary, db="gds", id=ids[0])[0]
    # el termino [ACCN] puede traer series relacionadas, no solo la pedida
    found = summary.get("Accession")
    if found != accession:
        raise ValueError(f"Se pidió {accession} pero NCBI devolvió {found}")
    return summary


def validate_bioproject(summary: dict[str, Any]) -> dict[str, bool]:
    """Busca el BioProject enlazado a la serie.

    Parameters
    ----------
    summary : dict
        Ficha devuelta por :func:`get_gse_summary`.

    Returns
    -------
    dict of {str: bool}
        Códigos PRJ encontrados como claves; vacío si no hay enlace.

    Warns
    -----
    UserWarning
        Si GEO no enlaza un BioProject o la respuesta no se puede leer.
    """
    links = entrez_query(
        Entrez.elink, dbfrom="gds", db="bioproject", id=summary["Id"]
    )
    ids = [
        link["Id"]
        for link_set in links
        for database in link_set.get("LinkSetDb", [])
        for link in database["Link"]
    ]
    # el BioProject es una comprobacion extra, no debe frenar el pipeline
    if not ids:
        warnings.warn(
            "GEO no enlaza un BioProject a esta serie vía Entrez; "
            "compruébalo a mano en la página de la serie.",
            UserWarning,
        )
        return {}
    result = entrez_query(
        Entrez.esummary, db="bioproject", id=",".join(ids)
    )
    try:
        docs = result["DocumentSummarySet"]["DocumentSummary"]
        return {doc["Project_Acc"]: True for doc in docs}
    except (KeyError, TypeError):
        warnings.warn(
            f"No se pudieron leer los códigos PRJ; se recibió: {result}",
            UserWarning,
        )
        return {}


def get_samples_table(summary: dict[str, Any]) -> pd.DataFrame:
    """Convierte la lista de muestras de la ficha en una tabla.

    Parameters
    ----------
    summary : dict
        Ficha devuelta por :func:`get_gse_summary`.

    Returns
    -------
    pandas.DataFrame
        Columnas ``gsm`` y ``title``, una fila por muestra.

    Raises
    ------
    ValueError
        Si la ficha no trae muestras o no coinciden con ``n_samples``.
    """
    try:
        rows = [
            {"gsm": sample["Accession"], "title": sample["Title"]}
            for sample in summary["Samples"]
        ]
        expected = int(summary["n_samples"])
    except KeyError as error:
        raise ValueError(
            f"La ficha no tiene el campo {error}. Solución: usa la salida "
            "de get_gse_summary()."
        ) from error
    # una lista incompleta cambiaria la seleccion de muestras despues
    if len(rows) != expected:
        raise ValueError(f"Se esperaban {expected} muestras y hay {len(rows)}")
    return pd.DataFrame(rows)


def list_supplementary_files(
    summary: dict[str, Any],
) -> tuple[str, list[str]]:
    """Lista los archivos de la carpeta ``suppl/`` de la serie.

    Parameters
    ----------
    summary : dict
        Ficha devuelta por :func:`get_gse_summary`.

    Returns
    -------
    suppl_url : str
        URL HTTPS de la carpeta ``suppl/``.
    files : list of str
        Nombres de los archivos disponibles.

    Raises
    ------
    FileNotFoundError
        Si la carpeta no existe en el servidor.
    ConnectionError
        Si no hay conexión con el servidor.
    """
    # urllib no lista carpetas FTP, pero el mismo servidor las sirve en HTTPS
    ftp = summary["FTPLink"].replace("ftp://", "https://")
    suppl_url = ftp.rstrip("/") + "/suppl/"
    try:
        with urllib.request.urlopen(suppl_url, timeout=60) as response:
            html = response.read().decode()
    except urllib.error.HTTPError as error:
        raise FileNotFoundError(
            f"No se pudo abrir {suppl_url} (HTTP {error.code}). Solución: "
            "revisa que la serie tenga archivos suplementarios."
        ) from error
    except urllib.error.URLError as error:
        raise ConnectionError(
            f"Sin conexión con {suppl_url}. Solución: revisa tu internet."
        ) from error
    links = re.findall(r'href="([^"?/][^"]*)"', html)
    # ':' o '/' indican enlaces externos o carpetas, no archivos de la serie
    files = [link for link in links if ":" not in link and "/" not in link]
    return suppl_url, files


def remote_exists(url: str) -> bool:
    """Comprueba si un archivo existe en el servidor sin descargarlo.

    Parameters
    ----------
    url : str
        Dirección completa del archivo.

    Returns
    -------
    bool
        True si el servidor responde sin error.

    Raises
    ------
    ConnectionError
        Si no hay conexión con el servidor.
    """
    # HEAD evita bajar el archivo completo solo para saber si existe
    request = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(request, timeout=60):
            return True
    except urllib.error.HTTPError:
        return False
    except urllib.error.URLError as error:
        raise ConnectionError(
            f"Sin conexión con {url}. Solución: revisa tu internet."
        ) from error


def is_gzip(path: str | Path) -> bool:
    """Comprueba que un archivo sea gzip por sus dos primeros bytes.

    Parameters
    ----------
    path : str or pathlib.Path
        Ruta del archivo.

    Returns
    -------
    bool
        True si empieza con la firma gzip ``1f 8b``.

    Raises
    ------
    FileNotFoundError
        Si el archivo no existe.
    """
    try:
        with open(path, "rb") as file:
            return file.read(2) == b"\x1f\x8b"
    except FileNotFoundError as error:
        raise FileNotFoundError(f"No existe: {path}") from error


def download_file(suppl_url: str, name: str,
                  dest_folder: str | Path) -> Path:
    """Descarga un archivo si existe en el servidor y aún no está local.

    Parameters
    ----------
    suppl_url : str
        URL de la carpeta ``suppl/`` (termina en ``/``).
    name : str
        Nombre del archivo a descargar.
    dest_folder : str or pathlib.Path
        Carpeta destino; se crea si no existe.

    Returns
    -------
    pathlib.Path
        Ruta local del archivo.

    Raises
    ------
    FileNotFoundError
        Si el archivo no existe en el servidor.
    ValueError
        Si la descarga queda vacía o no es un gzip válido.
    ConnectionError
        Si la conexión falla durante la descarga.
    """
    folder = Path(dest_folder)
    folder.mkdir(parents=True, exist_ok=True)
    destination = folder / name
    # reutilizar la copia local evita repetir una descarga de varios MB
    if destination.exists() and destination.stat().st_size > 0:
        print(f"Ya existe: {destination}")
        return destination
    url = suppl_url + name
    if not remote_exists(url):
        raise FileNotFoundError(
            f"No existe en el servidor: {url}. Solución: elige un nombre "
            "de list_supplementary_files()."
        )
    try:
        urllib.request.urlretrieve(url, destination)
    except (urllib.error.URLError, http.client.HTTPException) as error:
        # una descarga a medias pasaria por buena en la siguiente ejecucion
        destination.unlink(missing_ok=True)
        raise ConnectionError(
            f"La descarga de {name} se interrumpió. Solución: vuelve a "
            "ejecutar la celda."
        ) from error
    size = destination.stat().st_size
    if size == 0 or (name.endswith(".gz") and not is_gzip(destination)):
        destination.unlink()
        raise ValueError(
            f"{name} llegó vacío o dañado y se borró. Solución: vuelve a "
            "ejecutar la celda para descargarlo de nuevo."
        )
    print(f"Descargado: {destination} ({size / 1e6:.2f} MB)")
    return destination


def validate_count_matrix(path: str | Path) -> pd.DataFrame:
    """Lee la matriz y comprueba que tenga formato de conteos crudos.

    Parameters
    ----------
    path : str or pathlib.Path
        Ruta de la matriz separada por tabuladores (admite ``.gz``).

    Returns
    -------
    pandas.DataFrame
        Matriz genes × muestras.

    Raises
    ------
    FileNotFoundError
        Si el archivo no existe.
    ValueError
        Si no se puede leer, tiene menos de 2 columnas o columnas no
        numéricas.

    Warns
    -----
    UserWarning
        Si hay valores no enteros, NaN o IDs duplicados.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"No existe: {path}. Solución: ejecuta antes download_file()."
        )
    try:
        matrix = pd.read_csv(path, sep="\t", index_col=0)
    except (pd.errors.ParserError, UnicodeDecodeError, EOFError) as error:
        raise ValueError(
            f"No se pudo leer {path.name}: {error}. Solución: borra el "
            "archivo y descárgalo de nuevo."
        ) from error
    # con una sola columna, el separador no era tabulador
    if matrix.shape[1] < 2:
        raise ValueError("La matriz tiene <2 columnas: revisa el separador")
    numeric = matrix.select_dtypes("number")
    if numeric.shape[1] != matrix.shape[1]:
        raise ValueError("Hay columnas no numéricas en la matriz")
    # DESeq2 exige conteos crudos; decimales indicarian datos normalizados
    all_integers = bool(((numeric.fillna(0) % 1) == 0).all().all())
    n_nan = int(matrix.isna().sum().sum())
    n_duplicated = int(matrix.index.duplicated().sum())
    print(f"Genes x muestras: {matrix.shape}")
    print(f"¿Todos enteros? {all_integers}")
    print(f"NaN: {n_nan} | IDs duplicados: {n_duplicated}")
    if not all_integers:
        warnings.warn(
            "La matriz no parece de conteos crudos; revísala antes de "
            "seguir.",
            UserWarning,
        )
    if n_nan or n_duplicated:
        warnings.warn(
            "Hay NaN o IDs duplicados; build_count_table() los rechazará.",
            UserWarning,
        )
    return matrix
