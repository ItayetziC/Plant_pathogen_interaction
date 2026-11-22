# Análisis transcriptómico diferencial de Pseudomonas syringae durante la infección de Arabidopsis thaliana

## Se ubica en [Plant_pathogen_interaction](https://github.com/ItayetziC/Plant_pathogen_interaction)

## Funcionalidades
-  Descarga y limpieza de matriz de conteo de expresión de genes condiciones de infección de _Pseudomonas syringae_ (_in planta_) frente a controles de crecimiento no infectivo (_in vitro_). 
- Descarga y valida los registros de GenBank para la anotación funcional de la bacteria, a partir del genoma secuenciado por Buell et al. (2003).
- Control de calidad de datos por medio de un PCA de las muestras obtenidas.
- Análisis de expresión diferencial con InMoose: normalización por mediana de razones, estimación de dispersión con encogimiento bayesiano, modelo binomial negativo y prueba de Wald con corrección de Benjamini-Hochberg.
- Validación del análisis contra genes control de comportamiento conocido.
- Visualización de resultados: volcano plot, mapa de calor de los genes más expresados y figura de candidatos sin función asignada.



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
- Las bibliotecas listadas en requirements.txt: biopython, inmoose, pandas, numpy, scipy, statsmodels, matplotlib y seaborn.
- Conexión a internet para consultar NCBI Entrez, descargar los datos de GEO y consultar la base de datos KEGG.
- Una dirección de correo electrónico válida: NCBI lo necesita para identificar las consultas a Entrez. Se configura en la celda 1 del notebook 01_descarga.ipynb (Entrez.email)
  antes de ejecutarlo.
- Menos de 1 GB de memoria RAM y unos 5 MB de espacio en disco para los datos descargados.
Nota adicional: La instalación de InMoose se distribuye únicamente como código fuente y contiene extensiones en C++, por lo que su instalación requiere un compilador.
En Windows hay que instalar previamente las Microsoft C++ Build Tools con el componente "Desarrollo de escritorio con C++"; en Linux y macOS el compilador suele estar disponible.


## Datos
El proyecto utiliza dos fuentes públicas del NCBI: 
- **Matriz de conteos**: serie [GSE103441](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103441) de GEO (Nobori et al., 2018), con 114 muestras de RNA-seq bacteriano _in planta_ y en medio de cultivo. 
- **Genoma de referencia**: registros `NC_004578.1` (cromosoma), `NC_004632.1` (plásmido pDC3000A) y `NC_004633.1` (plásmido pDC3000B) de NCBI Nucleotide. 

Los datos no se versionan en este repositorio: se descargan automáticamente al ejecutar el notebook. Antes de hacerlo, configura la dirección de correo en la celda 1 (`Entrez.email`), ya que NCBI la exige para identificar las consultas. 
Los archivos descargados se guardan en `data/raw/` y los archivos derivados en `data/processed/`. Ambas carpetas se crean automáticamente si no existen. La procedencia, el diseño experimental y los criterios de selección de muestras se documentan en el [reporte del proyecto](docs/project_report.md).


La procedencia y características detalladas se describen en el
[reporte del proyecto](docs/reporte-proyecto.md).

## Uso

El flujo se ejecuta desde cinco notebooks numerados que importan las funciones de src/. 
Deben ejecutarse en orden, porque cada uno utiliza los archivos generados del anterior. 
Desde la raíz del repositorio, se muestra un ejemplo:

```bash
jupyter notebook notebooks/01_descarga.ipynb
```
El comando generará la descarga de datos de Geo y registros de GenBank seleccionados (ver _docs/project_report_) y limpia muestras depositados en _data/raw_ y _data/processed_.  

|Notebook |	Qué hace |
|:--|:--|
| 01_descarga.ipynb	| Descarga los datos de GEO y los registros GenBank, selecciona y limpia las muestras, construye la anotación funcional y evalúa la calidad con un PCA |
| 02_DEGs.ipynb	| Normaliza, estima la dispersión, ajusta el modelo binomial negativo, aplica la prueba de Wald y guarda la tabla de resultados |
| 03_graficas.ipynb	| Genera el volcano plot, el mapa de calor y la figura de candidatos sin función asignada |
| 04_KEGG.ipynb |	Obtiene las asociaciones gen-ruta de KEGG y ejecuta el análisis de sobrerrepresentación |
| 05_analisis_funcional.ipynb |	Clasifica los genes del sistema de secreción tipo III e identifica los de función desconocida |

### Scripts usados
Las funciones están en src/ y los notebooks solo las llaman. 
Cada función está documentada con su propósito, sus parámetros de entrada y su valor de retorno, y valida sus entradas lanzando excepciones con un mensaje que indica cómo corregir 
el problema.

| Módulo | Responsabilidad |
|:--|:--|
| geo_download.py |	Consulta GEO con Entrez y descarga los archivos de la serie |
| data_prep.py	| Valida la matriz, selecciona las muestras del contraste, limpia y filtra |
| annotation.py	| Descarga los registros GenBank y une la anotación a la matriz de conteos |
| quality_control.py	| Métricas de profundidad por muestra y análisis de componentes principales |
| differential_expression.py	| Construcción del modelo, normalización, dispersión, prueba de Wald y clasificación |
| kegg_analysis.py	Consulta de rutas KEGG y análisis de sobrerrepresentación |
| functional_analysis.py	| Clasificación funcional de los genes diferencialmente expresados |
| visualization.py	| Todas las figuras del proyecto |


## Reproducción de resultados

Antes de ejecutar, modifica Entrez.email en la celda 1 de 01_descarga.ipynb. Después ejecuta los cinco notebooks en orden: 

```bash
git clone https://github.com/ItayetziC/Plant_pathogen_interaction.git
cd Plant_pathogen_interaction
pip install -r requirements.txt
jupyter notebook notebooks/
```

Al terminar el notebook se habrán generado los siguientes archivos:

### Archivos generados: 

*Datos (ubicados en _/data_)*

| Archivo	| Contenido |
|:--|:--|
| raw/GSE103441_count.txt.gz	| Matriz de conteos original descargada de GEO |
| raw/NC_004578.1.gb y los dos plásmidos	| Registros GenBank del genoma de DC3000 |
| processed/GEO_samples.csv	| Las 114 muestras de la serie con su código GSM |
| processed/counts_filtered.csv	| Matriz de conteos crudos filtrada, restringida a las muestras del contraste |
| processed/sample_metadata.csv |	Condición experimental y código GSM de cada muestra seleccionada |
| processed/genbank_annotation_DC3000.csv	| Anotación extraída de los registros GenBank |
| processed/final_annotation.csv	| Anotación funcional por gen, combinando GenBank y la matriz |
| data/processed/genbank_annotation_DC3000	.csv	|	Anotación funcional para todo el genoma de GenBank |
| results/figures/library_metrics.png	|	Gráficas de Profundidad (millones de lecturas) y Mediana de conteos por gen |
La ejecución completa tarda unos pocos minutos y requiere conexión a internet.

*Tablas de resultados (Ubicados en _results/tables/_)*
| Archivo |	Contenido |
|:--|:--|
| factores_de_tamano.csv	| Factor de tamaño estimado para cada muestra |
| deg_results_completa.csv	| Resultado del modelo para todos los genes evaluados, con su anotación |
| deg_results_significativos.csv	| Subconjunto de genes diferencialmente expresados |
| validacion_controles.csv	| Comportamiento de los genes control frente a lo esperado |
| kegg_pathways.tsv y kegg_gene_pathway.tsv	| Catálogo de rutas KEGG y su asociación con cada gen |
| ora_significant_pathways.csv	| Rutas enriquecidas entre los genes diferencialmente expresados |
| degs_funcion_desconocida.csv	| Genes diferenciales anotados como proteína hipotética |

*Figuras (Ubicados en _results/figures/_)*
| Archivo	| Contenido |
|:--|:--|
| library_metrics.png	| Profundidad y cobertura por muestra |
| pca_samples.png	| Análisis de componentes principales de las muestras |
| volcano_plot.png	| Magnitud y significancia del cambio para todos los genes |
| heatmap_top40.png	| Los 20 genes más inducidos y los 20 más reprimidos, por muestra |
| candidatos_hipoteticos.png	| Genes sin función asignada con cambio de expresión relevante |
| kegg_enrichment_dotplot.png	| Rutas KEGG enriquecidas |

## Resultados más relevantes 
a) A partir del Análisis de Expresión Diferencial con _Inmoose_:
1. De 5394 genes evaluados, 1458 resultaron diferencialmente expresados (FDR < 0.05 y |log2FoldChange| > 1): 830 inducidos y 628 reprimidos durante la infección.
2. Los cuatro controles positivos de virulencia se indujeron como predice la literatura: pvdS (+6.86), avrPto1 (+6.05), cmaA (+5.68) y hrpL (+4.90).
3. Entre los genes anotados como proteína hipotética, 122 presentan expresión diferencial. 

b) A partir del análisis de _ORA_:

La interpretación completa está en el reporte del proyecto, además se puede visualizar los resultados más relevantes en una presentación de [Canva](https://canva.link/6baxis1zo5ls3rw).

## Documentación

Para más información consulta:
- [Reporte del proyecto](docs/project_report.md) -> contiene el contexto biológico, metodología más detallada y muestra resultados de forma más extensa
- [Información para citar el software](CITATION.cff) -> en caso de citar este proyecto, por favor consúltelo 

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
