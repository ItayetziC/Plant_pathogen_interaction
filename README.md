# Análisis transcriptómico diferencial de Pseudomonas syringae durante la infección de Arabidopsis thaliana

# Se ubica en [Plant_pathogen_interaction](https://github.com/ItayetziC/Plant_pathogen_interaction)

## Funcionalidades

- [Funcionalidad disponible 1]
- [Funcionalidad disponible 2]
- [Funcionalidad disponible 3]

<!--
AYUDA:
Describan solamente capacidades que ya funcionan. Comiencen cada elemento con un
verbo. Las funciones planeadas deben registrarse como issues y/o en GitHub Projects.

EJEMPLO:
- Descarga genomas mediante identificadores de NCBI.
- Lee archivos FASTA y GFF3.
- Genera una tabla de presencia y ausencia de genes.
-->

## Estructura del repositorio

```text
Plant_pathogen_interaction/
├── data/           # Datos de ejemplo o archivos pequeños
├── docs/           # Reporte y documentación
├── notebooks/      # Análisis exploratorios
├── results/        # Tablas y figuras generadas
├── src/            # Código fuente
├── tests/          # Pruebas
├── CITATION.cff    # Información para citar el software
├── LICENSE         # Licencia
└── README.md       # Introducción y guía rápida de uso
```

<!--
AYUDA:
Modifiquen el árbol para que represente su repositorio real y eliminen carpetas
que no utilicen. Expliquen solamente los elementos principales; no enumeren cada
archivo. Eviten almacenar datos grandes, usuarios/contraseñas o información sensible.
-->

## Requisitos

- [Lenguaje y versión]
- [Herramienta o biblioteca indispensable]
- [Recurso computacional o condición de acceso]

<!--
AYUDA:
Indiquen lo necesario antes de instalar o ejecutar el proyecto: versión de Python,
sistema operativo si es relevante, memoria, almacenamiento o acceso a servicios.

EJEMPLO:
- Python 3.12 o posterior.
- Git.
- 4 GB de memoria RAM.
- Conexión a internet para obtener datos.

No registren usuarios, contraseñas, tokens ni llaves privadas-->


## Datos

[Expliquen cómo obtener los datos y dónde colocarlos.]

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

[Expliquen qué hace el siguiente comando y qué entradas necesita.]

```bash
python src/main.py [argumentos]
```

[Indiquen dónde se guardan los resultados.]

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

Ejecuten los siguientes pasos en el orden indicado:

```bash
python scripts/download_data.py
python src/run_analysis.py
python src/create_figures.py
```

Los resultados esperados se generarán en `results/`.

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

- [Reporte del proyecto](docs/reporte-proyecto.md)
- [Información para citar el software](CITATION.cff)

<!--
AYUDA:
Incluyan únicamente documentos que existan y comprueben sus enlaces. Agreguen
otros documentos sólo si evitan que el README sea demasiado extenso.

EJEMPLO:
- El reporte contiene el problema, las preguntas, la metodología y los resultados.
-->

## Equipo

- Andrómeda Coral Gutiérrez — PENDIENTE
- Itayetzi Corona Aquino — PENDIENTE
- Leilani Cruz Ramírez — PENDIENTE

<!--
AYUDA:
Identifiquen a las tres personas. Describan brevemente sus contribuciones generales.
Las actividades específicas y revisiones se consultan en los issues, Pull Requests,
GitHub Projects y el historial de Git.

EJEMPLO:
- Ana Pérez — procesamiento de datos y pruebas.
- Luis López — análisis y visualización.
- María García — documentación e integración.
-->

## Citación

Si utilizas este software, consulta [CITATION.cff](CITATION.cff) o la opción
**Cite this repository** de GitHub.

<!--
AYUDA:
Mantengan CITATION.cff actualizado con autores, título y versión. Si el proyecto
obtiene un DOI, añadan aquí la referencia recomendada.

EJEMPLO:
La forma recomendada de citar la versión v1.0.0 se encuentra en CITATION.cff.
-->

## Licencia

El código se distribuye bajo la licencia MIT. Consulta el [LICENSE](https://github.com/ItayetziC/Plant_pathogen_interaction/blob/main/LICENSE) para conocer
los términos de uso.

<!--
AYUDA:
Especifiquen la licencia del código y asegúrense de incluir el archivo LICENSE.
Las condiciones de uso de los datos pueden ser distintas y deben indicarse en la
sección Datos y en el reporte.

EJEMPLO:
El código se distribuye bajo la licencia MIT. Consulta LICENSE.
-->

## Agradecimientos

[Incluyan reconocimientos institucionales o académicos]

<!--
AYUDA:
Esta sección es opcional. Reconozcan laboratorios, docentes, instituciones,
proyectos o financiamientos que apoyaron el trabajo. No incluyan como autores a
personas que sólo deban aparecer en los agradecimientos.

EJEMPLO:
Proyecto desarrollado en la Licenciatura en Ciencias Genómicas, UNAM, como parte
de la asignatura [nombre].
-->

---

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
