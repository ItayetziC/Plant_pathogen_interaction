# Análisis transcriptómico diferencial de Pseudomonas syringae durante la infección de Arabidopsis thaliana

## Se ubica en [Plant_pathogen_interaction](https://github.com/ItayetziC/Plant_pathogen_interaction)

## Funcionalidades
-  Descarga y limpieza de matriz de conteo de expresión de genes condiciones de infección de _Pseudomonas syringae_ (_in planta_) frente a controles de crecimiento no infectivo (_in vitro_). 
- Descarga y valida los registros de GenBank para la anotación funcional de la bacteria, a partir del genoma secuenciado por Buell et al. (2003).
- Control de calidad de datos por medio de un PCA de las muestras obtenidas.


## Estructura del repositorio

```text
Plant_pathogen_interaction/
├── data/            # Datos descargados (raw) y derivados (processed); no versionados
├── docs/            # Reporte y documentación
├── notebooks/       # Análisis exploratorios
├── results/         # Tablas y figuras generadas
├── src/             # Código fuente
├── tests/           # Pruebas
├── CITATION.cff     # Información para citar el software
├── LICENSE          # Licencia
└── README.md        # Introducción y guía rápida de uso
└── requirements.txt # Lista de librerías requeridas

```



## Requisitos

- Python 3.13 o posterior.
- Git.
- Jupyter Notebook o JupyterLab para ejecutar `notebooks/01_descarga.ipynb`.
- Bibliotecas de Python: `biopython`, `pandas`, `numpy`, `matplotlib`.
- Conexión a internet para consultar NCBI Entrez y descargar los datos de GEO.
- Una dirección de correo electrónico válida: NCBI la exige para identificar
  las consultas a Entrez. Debe configurarse en la celda 1 del notebook
  (`Entrez.email`) antes de ejecutarlo.
- Menos de 1 GB de memoria RAM y unos 5 MB de espacio en disco para los datos
  descargados.


## Datos
El proyecto utiliza dos fuentes públicas del NCBI: 
- **Matriz de conteos**: serie [GSE103441](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103441) de GEO (Nobori et al., 2018), con 114 muestras de RNA-seq bacteriano _in planta_ y en medio de cultivo. 
- **Genoma de referencia**: registros `NC_004578.1` (cromosoma), `NC_004632.1` (plásmido pDC3000A) y `NC_004633.1` (plásmido pDC3000B) de NCBI Nucleotide. 

Los datos no se versionan en este repositorio: se descargan automáticamente al ejecutar el notebook. Antes de hacerlo, configura la dirección de correo en la celda 1 (`Entrez.email`), ya que NCBI la exige para identificar las consultas. 
Los archivos descargados se guardan en `data/raw/` y los archivos derivados en `data/processed/`. Ambas carpetas se crean automáticamente si no existen. La procedencia, el diseño experimental y los criterios de selección de muestras se documentan en el [reporte del proyecto](docs/project_report.md).

<!--
AYUDA:
Indiquen la fuente, el comando o enlace de descarga y la carpeta de destino. No
repitan aquí la descripción completa de muestras, formatos, variables o criterios
de selección; enlacen el reporte. Si los datos no pueden publicarse, expliquen los
requisitos y el procedimiento autorizado para acceder a ellos.

EJEMPLO:
Los datos proceden de NCBI RefSeq. Para descargarlos, ejecuta:

    python scripts/download_data.py

Los archivos se guardarán en data/raw/. Los identificadores, versiones y criterios
de selección se documentan en docs/reporte-proyecto.md.
-->

La procedencia y características detalladas se describen en el
[reporte del proyecto](docs/reporte-proyecto.md).

## Uso

El pipeline se ejecuta desde un notebook de Jupyter que importa las funciones
de `src/`. Desde la raíz del repositorio:

```bash
jupyter notebook notebooks/01_descarga.ipynb
```

Ejecuten las celdas en orden. El notebook realiza tres tareas encadenadas:
consulta y descarga los datos de GEO, selecciona y limpia las muestras de
interés, y construye la anotación funcional del genoma a partir de los
registros GenBank.

No requiere argumentos de línea de comandos. El único valor que debe
modificarse antes de ejecutar es `Entrez.email` en la celda 1.

Los archivos generados se escriben en `data/processed/`  y `data/raw/` y las figuras en
`results/figures/`.


<!--
AYUDA:
Incluyan al menos un ejemplo mínimo que pueda copiarse y ejecutarse. Sustituyan
los corchetes por valores reales. Expliquen entradas, opciones importantes y
archivos de salida sin describir toda la metodología.

EJEMPLO:
Para analizar los identificadores incluidos en data/accessions.txt:

    python src/main.py --input data/accessions.txt --output results/

El comando generará results/gene_matrix.csv y results/heatmap.png.
-->


## Reproducción de resultados

Desde un entorno limpio, ejecuten los siguientes pasos en orden:

```bash
git clone https://github.com/ItayetziC/Plant_pathogen_interaction.git
cd Plant_pathogen_interaction
pip install -r requirements.txt
jupyter notebook notebooks/01_descarga.ipynb
```

Antes de ejecutar las celdas, modifiquen `Entrez.email` en la celda 1.

Al terminar el notebook se habrán generado los siguientes archivos:

| Archivo | Contenido |
|:--|:--|
| `data/raw/GSE103441_count.txt.gz` | Matriz de conteos original descargada de GEO |
| `data/raw/NC_004578.1.gb` y los dos plásmidos | Registros GenBank del genoma de DC3000 |
| `data/processed/GEO_samples.csv` | Las 114 muestras de la serie con su código GSM |
| `data/processed/counts_filtered.csv` | Matriz de conteos crudos filtrada, restringida a las muestras del contraste |
| `data/processed/sample_metadata.csv` | Condición experimental y código GSM de cada muestra seleccionada |
| `data/processed/final_annotation.csv` | Anotación funcional por gen, combinando GenBank y la matriz |
| `results/figures/pca_samples.png` | Análisis de componentes principales de las muestras |

La ejecución completa tarda unos pocos minutos y requiere conexión a internet.
<!--
AYUDA:
Proporcionen la ruta más corta para regenerar el resultado principal desde los
datos originales. Los comandos deben indicar el orden correcto. Si el proceso
requiere parámetros o archivos de configuración, indíquenlos.
Las explicaciones científicas e interpretación de resultados pertenecen al reporte.

EJEMPLO:
Después de ejecutar los tres comandos se crearán la tabla comparativa y las
figuras utilizadas en el reporte. Sus nombres esperados deben indicarse aquí.
-->


## Documentación

Para más información consulta:
- [Reporte del proyecto](docs/project_report.md)
- [Información para citar el software](CITATION.cff)

<!--
AYUDA:
Incluyan únicamente documentos que existan y comprueben sus enlaces. Agreguen
otros documentos sólo si evitan que el README sea demasiado extenso.

EJEMPLO:
- El reporte contiene el problema, las preguntas, la metodología y los resultados.
-->

## Equipo

- Andrómeda Coral Gutiérrez — Descarga y procesamiento de datos e integración a Git
- Itayetzi Corona Aquino — Análisis bioinformático y visualización
- Leilani Cruz Ramírez — Documentación e interpretación biológica de los datos


## Citación

Si utilizas este software, consulta [CITATION.cff](CITATION.cff) o la opción
**Cite this repository** de GitHub.


## Licencia

El código se distribuye bajo la licencia MIT. Consulta el [LICENSE](https://github.com/ItayetziC/Plant_pathogen_interaction/blob/main/LICENSE) para conocer
los términos de uso.


## Agradecimientos

Proyecto desarrollado en la Licenciatura en Ciencias Genómicas, UNAM, como parte de la asignatura de Biopython.
Se le agradece a las docentes; Dra. Verónica Jiménez, Dra. Leticia Vega y a la maestra Shirley Alquicira por su apoyo. Además se le agradece a Pablo Salazar por su asesoría durante la realización del proyecto. 


<!--
LISTA DE COMPROBACIÓN ANTES DE ENTREGAR:
- [ ] Se eliminaron o sustituyeron todos los textos entre corchetes.
- [ ] Los enlaces funcionan desde GitHub.
- [ ] Los comandos fueron probados en un ambiente limpio.
- [ ] El README permite instalar y ejecutar un ejemplo.
- [ ] No se incluyeron credenciales ni información sensible.
- [ ] El reporte contiene la explicación científica y evita duplicar el README.
- [ ] La versión, la licencia y la información de citación están actualizadas.
-->
