# Análisis transcriptómico diferencial de *Pseudomonas syringae* durante la infección de *Arabidopsis thaliana*

<!-- AYUDA: Escriban un título breve, específico y descriptivo. Puede ser provisional.
EJEMPLO: Comparación de genes de resistencia antimicrobiana en genomas de
Escherichia coli. -->

## Información general

| Dato | Información |
|:--|:--|
| Integrante 1 | Itayetzi Corona Aquino itayetzi@lcg.unam.mx |
| Integrante 2 | Leilani Cruz Ramírez leilanic@lcg.unam.mx |
| Integrante 3 | Andrómeda Coral Gutiérrez andcoral@lcg.unam.mx |
| Fecha de creación | 25/08/2026 |
| Última actualización | 23/11/2026 |
| Repositorio | https://github.com/ItayetziC/Plant_pathogen_interaction |

<!-- AYUDA: Mantengan actualizados el estado, la fecha y la versión. La versión
debe coincidir con una etiqueta de Git cuando exista una entrega identificable. -->

## Resumen del proyecto
Análisis transcriptómico diferencial en Python de _Pseudomonas syringae_ pv. tomato DC3000 durante la infeccion de _Arabidopsis thaliana_. 
Integra datos públicos de RNA-seq (GEO GSE103441), análisis estadístico de expresión diferencial mediante InMoose, anotación funcional (RefSeq GCF_000007805.1) 
y enriquecimiento de rutas biológicas con KEGG, así como herramientas de visualización para explorar los mecanismos moleculares asociados con la interacción 
planta-patógeno.

## 1. Contexto y antecedentes

<!-- AYUDA: Presenten la información necesaria para comprender el proyecto.
Definan conceptos biológicos y computacionales, describan qué se conoce y citen
trabajos, datos o herramientas relacionados.
PREGUNTAS: ¿Cuál es el fenómeno de interés? ¿Qué debe conocer quien lea el
reporte? ¿Qué métodos o herramientas se han utilizado anteriormente?
EJEMPLO: La resistencia antimicrobiana es un problema de salud pública. Aunque
existen bases especializadas, comparar varios genomas requiere integrar datos
procedentes de distintos archivos. -->

_Pseudomonas syringae_ pv._tomato_ DC3000 es una bacteria Gram-negativa patógena que infecta 
tanto al tomate como a la planta modelo _Arabidopsis thaliana_, lo que la ha convertido en
uno de los sistemas de referencia más estudiados para descifrar los mecanismos moleculares de la 
patogénesis bacteriana y la inmunidad vegetal [1,2].

Para infectar a la planta, la bacteria ingresa a la hoja por los estomas, suprime la 
inmunidad vegetal para poder multiplicarse en el apoplasto y convertirlo en un medio acuoso 
(inundación del espacio intercelular) para absorber nutrientes [2,3]. 

El control transcripcional del proceso de infección ha sido ampliamente estudiado para 
comprender la red de genes virulentos, descritos gracias a la secuenciación completa del 
genoma de DC3000 que también reveló un cromosoma de 6.4 Mb y dos plásmidos que codifican 
5,763 marcos de lectura abierta, entre ellos 298 genes de virulencia establecidos y 
putativos [4]. 

El programa transcripcional que la bacteria despliega realmente dentro del hospedero permanece
menos caracterizado. Esto se debe a que, en tejido vegetal infectado, el ARN bacteriano representa una fracción 
mínima del ARN total, el cual está dominado por los transcritos de la planta. Sin embargo, en 2018 Nobori _et al_. 
[5] desarrollaron métodos de enriquecimiento que permitieron superar este obstáculo y 
generaron el transcriptoma de DC3000 durante la infección de _A. thaliana_, depositado públicamente 
como GSE103442 y GSE103441.

En este proyecto se realizará un análisis de expresión génica diferencial utilizando el 
conjunto de datasets de Buell _et al._, (2003) y Nobori _et al._, (2018) para comparar 
estadísticamente la cantidad de ARN producido por los genes (su nivel de expresión), en dos 
condiciones distintas: durante la infección a _A. thaliana_ y en condiciones no infectivas 
(crecimiento en medio de cultivo) [4,5]. 

Aproximadamente 658 genes estaban anotados como proteínas hipotéticas, lo que sugiere que una cantidad 
considerable de genes activos durante la infección carece de función asignada [5]. De modo
que se podrían identificar cambios de expresión cuyo significado biológico permanece indeterminado.

El análisis se implementará en Python con inmoose [6], y la anotación funcional de los genes 
se construirá con Biopython a partir de los registros GenBank del genoma (AE016853–AE016855). 
El resultado principal será una tabla de genes diferencialmente expresados con su magnitud de 
cambio [7] y su significancia estadística corregida, representada visualmente mediante un heatmap 
y un volcano plot [8].


## 2. Planteamiento del problema

<!-- AYUDA: Describan la dificultad, necesidad o vacío de conocimiento que
desean atender, a quién afecta y qué sucede si no se resuelve. No confundan el
problema con la herramienta ni con la solución.
EJEMPLO: La identificación manual de genes de resistencia en varios genomas es
lenta, propensa a errores y difícil de reproducir. -->


_Pseudomonas syringae_ pv. tomato DC3000 causa la mancha bacteriana del tomate, una enfermedad de 
importancia económica mundial para la cual aún no existen medidas de control efectivas [4]. 
Diseñarlas requiere saber qué genes utiliza la bacteria durante la infección real; sin embargo, buena parte 
de lo que se conoce sobre la regulación de su virulencia proviene de experimentos en medios de cultivo, 
los cuales no reproducen las condiciones del apoplasto. Además solo el 61% de sus 5763 genes tiene una función 
asignada, por lo que son necesarios más anáisis para comprender su mecanismo completo de patogenicidad.


## 3. Justificación

<!-- AYUDA: Expliquen por qué vale la pena realizar el proyecto, cuál es su
relevancia biológica, científica, técnica o social y quién podría beneficiarse.
EJEMPLO: Un flujo automatizado reducirá errores y permitirá repetir el análisis
con los mismos datos, parámetros y versiones del software. -->

_Pseudomonas_ es la causante de más de 50 tipos de enfermedades en diferentes cultivos que afectan 
tanto económicamente como socialmente en todo el mundo. 
Además, aún no existe una estrategia biotecnológica para controlar esta bacteria ni se conocen todos 
los genes involucrados en su patogenicidad que podrían ser clave para su desarrollo. 
Aunque los datos transcriptómicos necesarios para reducir este vacío ya son públicos, continúan 
subutilizados. Determinar cuáles cambian de forma significativa podría elucidar genes candidatos para priorizar en estudios funcionales, según su variación de su expresión al infectar a la planta 
modelo _A. thaliana_.

## 4. Objetivo general

<!-- AYUDA: Expresen el resultado global mediante un verbo en infinitivo. Debe
ser alcanzable durante el semestre.
EJEMPLO: Desarrollar un flujo reproducible en Python para identificar y comparar
genes de resistencia en un conjunto de genomas de E. coli. -->

Identificar y cuantificar la expresión de los genes de _Pseudomonas syringae_
pv.tomato DC3000 durante la infección de _Arabidopsis thaliana_, en comparación con condiciones
de crecimiento no infectivas, para determinar las funciones biológicas asociadas con la 
interacción bacteria-hospedero.


## 5. Preguntas de investigación

<!-- AYUDA: Formulen preguntas biológicas o computacionales que puedan
responderse con los datos y métodos disponibles. Indiquen qué evidencia sería
necesaria.
EJEMPLO: ¿Qué genes de resistencia aparecen en cada genoma? Evidencia: una tabla
de presencia y ausencia obtenida de las anotaciones. -->

### Pregunta 1

**Pregunta:** ¿Qué genes de _Pseudomonas syringae_ pv. _tomato_ DC3000 presentan cambios 
significativos en su expresión durante la infección de _Arabidopsis thaliana_ en comparación 
con condiciones de crecimiento no infectadas?  

**Evidencia necesaria:** Matriz de expresión génica y réplicas experimentales de ambas condiciones,
a partir de las cuales se obtenga una tabla de genes diferencialmente expresados con sus valores 
y significancia. Los resultados podrán visualizarse mediante un volcano plot y un heatmap.


### Pregunta 2

**Pregunta:** ¿Qué funciones biológicas y procesos asociados con la adaptación al hospedero y la 
patogénesis están representados entre los genes diferencialmente expresados durante la infección de 
_Pseudomonas syringae_?

**Evidencia necesaria:** Anotación funcional de los genes diferencialmente expresados y clasificación 
de estos según sus funciones, particularmente aquellas relacionadas con secreción y efectores, transporte 
y adquisición de nutrientes, regulación, adhesión, respuesta al estrés y otros mecanismos asociados con la virulencia.

### Pregunta 3

**Pregunta:** ¿Existen genes diferencialmente expresados durante la infección de _Pseudomonas syringae_
cuya función se encuentre poco caracterizada o sea desconocida y que, por su patrón de expresión, puedan 
considerarse candidatos para estudios posteriores de la interacción bacteria-hospedero?

**Evidencia necesaria:** Identificación de genes diferencialmente expresados con anotaciones 
de función desconocida, hipotética o conservada, considerando su magnitud de cambio 
de expresión y significancia estadística.


## 6. Alcance y limitaciones

<!-- AYUDA: Delimiten organismos, muestras, datos, análisis y resultado esperado. Indiquen
qué no se abordará y las restricciones de tiempo, cómputo, acceso o calidad.
EJEMPLO: Se analizarán como máximo 20 genomas completos de RefSeq. No se
utilizarán datos clínicos ni se realizará validación experimental. -->

### Incluye
- **Análisis bioinformático** Control de calidad de datos, cuantificación de la expresión génica
- **Muestras y datos:** 
Datos de RNA-seq: (Serie GEO `GSE103441` / BioProject `PRJNA401856`, Nobori et al., 2018) que abarcan condiciones de infección (_in planta_) frente a controles de crecimiento no infectivo (_in vitro_).
- **Genoma de referencia:** 
Ensamblado de RefSeq `GCF_000007805.1` para el alineamiento y anotación funcional de la bacteria realizado por Buell et al., (2003)

### No incluye

- **Validación experimental** (Ensayos de laboratorio (_wet-lab_), como RT-qPCR, mutagénesis ni infecciones _in vivo_).
- **Otras cepas u organismos:** No se evaluarán otras cepas de _P. syringae_ ni otros patógenos u 
hospederos alternativos.    
-   **Transcriptómica del hospedero:** No se caracterizará la respuesta transcripcional de 
*Arabidopsis thaliana*, enfocándose únicamente en el transcriptoma bacteriano.
-   **Secuenciación de nuevo (*de novo*):** No se generarán datos de secuenciación propios ni se 
realizará ensamblado _de novo_.

### Limitaciones conocidas

- **Origen de los datos:** Dependencia total de la calidad, profundidad de secuenciación y diseño  
experimental de los datos públicos depositados en NCBI por Nobori _et al,._ (2018) y Buell _et al,._ (2003).
	- **Número de réplicas:** De la bacteria se tienen 13 muestras _in planta_ vs. 2 muestras _in vitro_. Lo que baja el poder estadístico.
- **Falta de anotación certera:** A pesar de identificar cambios de expresión, no se podrá asegurar la función
- de genes putativos en el proceso de infección.

## 7. Propuesta de solución

<!-- AYUDA: Describan el producto o estrategia que podría resolver el problema.
Es una propuesta inicial y puede cambiar. Expliquen sus componentes, no sólo las
tecnologías.
EJEMPLO: Un programa modular recibirá identificadores, descargará archivos,
extraerá genes, almacenará resultados y generará visualizaciones. -->

1. Existe variación entre las muestras de la condición _in planta_ debido a la diferencia de la profundidad de las lecturas, es por ello que se implementará inmoose, donde internamente normaliza los datos antes de hacer el DEG.
2. La reducida cantidad de muestras control en el dataset: debido al análisis preliminar como control de calidad de los datos, se observó que estaban muy separadas las muestras
3. por lo que no se consideró necesario obtener más datos de otros datasets o otras condiciones que metieran más ruido al análisis, se optó por usar la dispersión y encogimiento bayesiano para mitigar (aunque no lo erradica) el problema de tener solo 2 muestras como control

### 7.1 Resultado o producto esperado

<!-- AYUDA: Indiquen el entregable concreto: programa, paquete, flujo de análisis,
base de datos, visualizaciones u otro producto.
EJEMPLO: Repositorio ejecutable con scripts, datos de prueba, documentación,
tabla comparativa y figuras regenerables. -->

1. Se espera un repositorio ejecutable por medio de scripts y notebooks que sea reproducible y tenga una documencación descriptiva.
   Además las figuras y resultados serán regenerables. Así como el uso de comprobaciones y aviso de errores explícito. 
2. Se descargarán los datos esperados, se realizará el análisis de expresión diferencial por medio de inmoose y para su visualización
   de usarán Volcano plot y Heatmaps
3. Para el enriquecimiento funcional con los genes con patrón de expresión diferencial se visualizará por medio de scatter plots.
4. Como parte de la obtención de genes sin anotación a candidatos virulentos se analizará su patrón de expresión y ubicación genómica,
   además de gráficos visuales que lo comprueben. 


## 8. Datos

### 8.1 Fuentes de datos

<!-- AYUDA: Incluyan institución, base de datos, URL, identificador, versión o
fecha de consulta y condiciones de uso. No todos los proyectos usarán NCBI.
EJEMPLO: NCBI RefSeq, GCF_000005845.2, consultado el dd/mm/aaaa. -->

| Fuente | Identificador o versión | URL | Fecha de consulta | Licencia o condiciones |
|:--|:--|:--|:--|:--|
| NCBI Assembly / RefSeq (NCBI / NIH) | Genoma de referencia *Pseudomonas syringae* pv. *tomato* DC3000	`GCF_000007805.1`  (NC_004578.1 , NC_004633.1 y NC_004632.1) |	[https://www.ncbi.nlm.nih.gov/assembly/29388](https://www.ncbi.nlm.nih.gov/assembly/29388)	|	30/08/2026	|	Dominio Público. Uso libre para investigación académica y comercial.
| NCBI GEO (Gene Expression Omnibus)	|	Datos de RNA-seq de *P. syringae* DC3000 (*in planta* vs. *in vitro*) Serie GEO: `GSE103441` (BioProject: `PRJNA401856`)	|	https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103441	|	30/08/2026	|	Datos de acceso abierto. Reutilización libre citando la publicación de Nobori _et al._ (2018).

### 8.2 Características de los datos

<!-- AYUDA: Describan organismos, muestras, variables, formatos, versiones, tamaño y otros
atributos necesarios para interpretar los datos.
EJEMPLO: Archivos FASTA y GFF3 de 20 genomas completos de E. coli. -->

**1. Descripción del genoma de referencia de _Pseudomonas syringae_**
-   **Organismo:** _Pseudomonas syringae_ pv. _tomato_ (cepa DC3000).

-   **Muestras / Tipo de datos:** Genoma bacteriano completo de referencia y su anotación funcional. 
-   **Variables / Atributos:** Secuencia completa de nucleótidos del cromosoma bacteriano 
(~6.39 Mb, ~5,600 genes codificantes de proteínas) y sus plásmidos asociados (pDC3000A y pDC3000B). 
Incluye coordenadas genómicas de exones, CDS, rRNAs, tRNAs y nombres de genes.  
-   **Formatos:** Archivo de secuencia en formato **FASTA** (`.fna`) y archivo de anotación genómica 
en formato **GFF3** (`.gff`) o **GTF**.
-   **Versión:** Ensamblado de RefSeq `GCF_000007805.1` (Secuencia de acceso al cromosoma principal: 
`NC_004578.1`; y Secuencias de acceso a los plásmidos: `NC_004633.1` y `NC_004632.1`).
-   **Tamaño:** ~2.0 MB comprimido (`.gz`).

**2.  Descripción de los datos transcriptómicos de *Pseudomonas syringae***

-   **Organismo:** *Pseudomonas syringae* pv. *tomato* (cepa DC3000).
-   **Muestras:** Muestras de perfiles transcriptómicos bacterianos derivados del procesamiento de la matriz de conteos (`GSE103441_count.txt.gz`), correspondientes a dos condiciones experimentales:
    1.  *In planta:* Bacterias aisladas o recolectadas del apoplasto de *Arabidopsis thaliana* 
    durante el proceso de infección.
    2.  *In vitro:* Bacterias cultivadas en medio sintético de laboratorio (condición control/no infectiva).
-   **Variables / Atributos:** Niveles de abundancia de transcritos de ARN (expresión génica bacteriana) 
en respuesta al entorno del hospedero. Mide la cantidad de lecturas mapeadas por gen bacteriano para 
comparar la activación de factores de virulencia, efectores y metabolismo adaptativo.

-   **Formatos:**
    -   **Datos procesados:** Matriz de conteos (*read counts table*) en texto plano delimitado 
    por tabuladores (`.txt` o `.tsv` comprimido en `.gz`), con filas representando IDs de genes de 
    DC3000 y columnas representando las réplicas biológicas.
    -   **Datos brutos:** Archivos de secuenciación de alto rendimiento (*Illumina*) en formato 
    **FASTQ** comprimido (`.fastq.gz`) almacenados en SRA.
-   **Versión:** Serie GEO `GSE103441` (Estudio de Nobori et al., 2018).
-   **Tamaño:** ~870.4 Kb para la tabla de datos procesados/conteos.


### 8.3 Organización de los datos

<!-- AYUDA: Muestren la estructura prevista. No suban datos sensibles, tokens,
contraseñas ni archivos grandes. Usen .gitignore y documenten cómo obtener lo
que no se guarde en Git.
EJEMPLO: data/raw conserva originales y data/processed los derivados. -->

```text
Plant_pathogen_interaction/
├── data/            # Datos descargados (raw) y derivados (processed); no versionados
├── docs/            # Reporte y documentación
├── notebooks/       # Análisis exploratorios
├── results/         # Tablas y figuras generadas; no versionados 
├── src/             # Código fuente
├── tests/           # Pruebas
├── CITATION.cff     # Información para citar el software
├── LICENSE          # Licencia
└── README.md        # Introducción y guía rápida de uso
└── requirements.txt # Lista de librerías requeridas
```
Para obtener los archivos no versionados: 
Ejecutar `notebooks/01_descarga.ipynb` en orden, de la celda 0 a la 13 (con las dependencias descritas en el README en reproducción de resultados) . Las carpetas `data/raw/`, `data/processed/` y `results/` se crean automáticamente si no existen.

### 8.4 Diccionario o formato de los datos

<!-- AYUDA: Describan campos o columnas relevantes. Incluyan fragmentos pequeños
cuando ayuden a comprender el formato, pero no archivos completos.
EJEMPLO: En GFF3, seqid identifica la secuencia; type indica gene, CDS, etc. -->

| Archivo o conjunto | Campo/columna | Tipo | Descripción | Valores o unidades |
|:--|:--|:--|:--|:--|
| `GSE103441_count.txt.gz` (GEO) | Primera columna (ID de gen) | Texto | Locus tag antiguo unido a la descripción del producto con un guion | `PSPTO_0001-chromosomal replication initiation protein` |
| `GSE103441_count.txt.gz` | Columnas de muestra | Entero ≥ 0 | Conteos crudos por gen; el nombre indica condición y réplica | `Pto_col_A13` (planta), `Pto_KB_A81` (King's B); lecturas. 5,618 genes × 114 muestras |
| GenBank `.gb` (NC_004578.1, NC_004633.1, NC_004632.1) | `CDS` y su ubicación | Rasgo | Región codificante con coordenadas y hebra | `CDS 339..1874`; `complement(...)` = hebra − |
| GenBank `.gb` | `/old_locus_tag` | Texto | Locus tag antiguo; clave para unir con la matriz | `PSPTO0001` |
| GenBank `.gb` | `/product` | Texto | Producto proteico anotado | `chromosomal replication initiator protein DnaA` |
| `sample_metadata.csv` | `condition` | Categórico | Condición de cada muestra elegida | `plant` (n=13), `culture` (n=2) |
| `final_annotation.csv` | `locus_key` | Texto | Locus tag sin guion bajo; une conteos y anotación | `PSPTO0001`; 5,394 genes |
| `final_annotation.csv` | `final_product`, `source` | Texto / categórico | Producto final y su fuente (GenBank si existe; si no, la matriz) | `GenBank (RefSeq)` 5,152; `Count matrix` 242 |

## 9. Metodología

<!-- AYUDA: Esta sección evolucionará. Primero describan el plan y después
actualícenla con lo que realmente ejecutaron, incluidos parámetros y decisiones. -->

### 9.1 Etapas del análisis o desarrollo

<!-- AYUDA: Describan la secuencia desde la obtención de datos hasta la validación
de resultados. Relacionen cada etapa con una pregunta u objetivo.
EJEMPLO: descarga, validación, transformación, análisis, visualización y pruebas. -->


El análisis se organiza en seis etapas. Las primeras tres preparan y validan los datos. Las últimas tres producen las respuestas a las preguntas de investigación.

**Etapa 1. Entrada.** El análisis parte de dos fuentes independientes. La primera es la matriz de conteos crudos de expresión de genes condiciones de infección de _Pseudomonas syringae_ (_in planta_) frente a controles de crecimiento no infectivo (_in vitro_) presente en la serie GSE103441 de GEO, con 114 muestras de RNA-seq bacteriano.  La segunda es el genoma anotado de *P. syringae* pv. tomato DC3000, compuesto por el cromosoma y sus dos plásmidos.

**Etapa 2. Preparación.** De las 114 muestras se seleccionan 15. Son las de la cepa silvestre en las dos condiciones que define el objetivo: 13 recuperadas de *A. thaliana* Col-0 y 2 crecidas en medio King's B. Archivos de GenBank en una tabla que asocia cada identificador con el producto de su gen. 

**Etapa 3. Control de calidad.** Se verifica la profundidad y la mediana de conteos de cada muestra. Se calcula un PCA que, sin recibir las etiquetas de condición, separa los dos grupos en el primer componente con un 43.8 % de la varianza. 

**Etapa 4. Análisis.** Los conteos crudos de las 15 muestras se modelarán con InMoose. 

**Etapa 5. Integración.** La tabla de resultados del modelo se unirá con la anotación funcional mediante el identificador de cada gen, dando como resultado una sola tabla que contiene, por gen, su magnitud de cambio, su significancia y la función de su producto. 

**Etapa 6. Interpretación.** De la tabla integrada se derivan tres salidas, una por cada pregunta:
	El volcano plot y el heatmap visualizan la magnitud y la dirección de los cambios, y completan la respuesta a la **Pregunta 1** (descriptiva). 
	El análisis de enriquecimiento identifica qué funciones biológicas aparecen sobrerrepresentadas entre los genes diferencialmente expresados, y responde la **Pregunta 2** (funcional). 
	El filtrado de los genes anotados como proteína hipotética que resultan inducidos durante la infección responde la **Pregunta 3** (exploratoria), al priorizar candidatos para caracterización funcional posterior.

**Figura de la estrategia metodológica** Se divide en 3 fases la obtención de datos (Etapa 1-3), Análisis de expresión diferencial (Etapa 4) e 
Interpretación (Etapa 5 y 6) 
<img width="1110" height="628" alt="imagen" src="https://github.com/user-attachments/assets/3137d5cb-53ab-4494-a640-af15b075c6a5" />

### 9.2 Herramientas y tecnologías

<!-- AYUDA: Registren lenguajes, bibliotecas y programas con sus versiones y
propósito. No incluyan credenciales.
EJEMPLO: Python 3.x; Biopython para leer formatos biológicos; Seaborn para
visualización. -->

| Herramienta | Versión | Propósito |
|:--|:--|:--|
| Python | 3.14.7 | Para la ejecución de código |
| Jupyter Notebook | 7.6.3 | Dividir el código en partes |
| Biopython | 1.88 | Acceder a Entrez | 
| Pandas | 2.3.3 | Manipulación y análisis de datos tabulares |
| Matplotlib.pyplot | 3.10.6 | Crear y modificar gráficos |
| Numpy | 2.5.2 | Manejo de datos y cálculos numéricos |
| Scipy | 1.16.3 | Cálculos matemáticos avanzados sobre arreglos de NumPy |
| Inmoose | 0.9.1 | Análisis de expresión diferencial |
| statsmodels | 0.14.5 | Uso obligatorio por Inmoose |
| seaborn | 0.13.2 | Gráficos visuales y estadísticos |


### 9.3 Estrategia de validación

Para garantizar la confiabilidad de los resultados y la reproducibilidad del análisis, se implementó una estrategia de validación que contempla la integridad de los datos de entrada, el control de calidad de las muestras, la consistencia de los análisis estadísticos y la verificación de las anotaciones funcionales.

***Validación de datos y funcionamiento del software:*** Se incorporaron 48 excepciones y comprobaciones en los scripts del proyecto para detectar entradas inválidas, archivos incompletos, identificadores inexistentes e inconsistencias en las estructuras de datos. Estas validaciones permiten identificar errores durante la ejecución y evitar que se propaguen hacia las etapas posteriores del análisis.

***Validación de muestras y control de calidad:*** Se evaluó la profundidad de secuenciación y la distribución de los conteos por muestra para identificar diferencias técnicas que pudieran afectar las comparaciones entre condiciones experimentales. Adicionalmente, se realizó un análisis de componentes principales (PCA), mediante el cual se observó una separación entre las muestras in planta e in vitro. El primer componente principal explicó el 43.8 % de la variabilidad, lo que proporciona evidencia exploratoria de diferencias entre ambos grupos.
El PCA confirmó la separación esperada y descartó réplicas atípicas; los análisis preliminares determinaron que se puede comparar por profundidad de las lecturas entre muestras a pesar de la dispersión in planta.

***Validación estadística del análisis de expresión diferencial:*** Se verificó la correspondencia entre los identificadores génicos y las anotaciones funcionales utilizadas para interpretar los resultados de expresión diferencial.La integración con KEGG se realizó mediante tablas que relacionan los identificadores de las rutas biológicas con sus nombres y los genes asociados a cada una.
El análisis de enriquecimiento funcional permite identificar rutas sobrerrepresentadas entre los genes diferencialmente expresados. Su significancia se evalúa mediante valores de probabilidad ajustados por comparaciones múltiples, considerando el conjunto de genes de referencia utilizado.

***Validación del modelo mediante controles biológicos:*** Se definieron dos conjuntos de genes con comportamiento conocido para comprobar que el análisis recupera la señal esperada. Como controles positivos se emplearon reguladores y factores de virulencia que deben inducirse durante la infección: *hrpL* (PSPTO1404), regulador maestro del sistema de secreción tipo III; *avrPto1* (PSPTO4001), efector; *cmaA* (PSPTO4709), biosíntesis de coronatina; y *pvdS* (PSPTO2133), adquisición de hierro. Como controles negativos se emplearon dos genes de mantenimiento celular: *gyrA* (PSPTO1745) y *gap* (PSPTO1287). El resultado de esta verificación se reporta en la sección 11 y se guarda en `results/tables/validacion_controles.csv`.

## 10. Plan de trabajo


### 10.1 Distribución de responsabilidades

<!-- AYUDA: Definan responsabilidades iniciales sin aislar a cada integrante.
Toda contribución importante debe ser revisada mediante Pull Request por otra
persona.
EJEMPLO: Ana desarrolla la descarga y revisa el módulo de visualización. -->

| Integrante | Responsabilidad principal | Responsabilidad de revisión |
|:--|:--|:--|
|   Itayetzi Corona Aquino  |   Codificación   | Análisis bioinformático |
|   Leilani Cruz Ramírez    |   Contexto biológico   | Interpretación biológica de los resultados |
|   Andrómeda Coral Gutiérrez   |   Supervisión    | Preparación de datos. Mantenimiento de repositorio |

### 10.2 Riesgos y alternativas

| Riesgo | Probabilidad |	Impacto	|	Prevención o alternativa
|	--	|	--	|	--	|	--	|	
|	Pocas réplicas en la condición de control (13 _in planta_ frente a 2 en King's B), lo que reduce el poder estadístico	|	Alta	|	Alto	|	Usar InMoose, cuyo encogimiento bayesiano estima la dispersión de cada gen combinando su variabilidad observada con la tendencia del conjunto de genes. Esto rescató 1448 genes cuya estimación individual resultaba imposible con dos réplicas. Complementariamente, se evaluó ampliar el control con cepas mutantes crecidas en King's B y se documentó la decisión de conservar el diseño original.
|	Desbalance entre condiciones y variación en la profundidad de secuenciación entre muestras _in planta_	|	Media	|	Medio	|	Normalizar con InMoose antes del análisis diferencial. Revisar la profundidad y la mediana de conteos en el control de calidad y descartar réplicas atípicas con el PCA.	
|	Dependencia total de la calidad y el diseño de datos públicos (Nobori _et al._, 2018), sin posibilidad de repetir o ampliar el experimento	|	Medio	|	Alto	|	Aclarar la fuente de los datos. Verificar el control de calidad del dataset. Apegarse a las condiciones del estudio
|	Anotación funcional incompleta: genes hipotéticos o sin función asignada, y desajustes al unir la matriz de conteos con GenBank	|	Alta	|	Medio	|	Unir por el identificador _locus tag_ antiguo, documentar los genes sin correspondencia y conservar la fuente de cada anotación. Presentar los genes hipotéticos como candidatos, sin asegurar su función.
|	Dificultades técnicas con la reproducibilidad (versiones de Python, librerías, o cambios en el formato de los archivos descargados)	|	Media	|	Medio	|		Fijar versiones en `project_report.md`, documentar la descarga en el notebook `01_descarga` y validar entradas con las excepciones ya implementadas.

## 11. Resultados

<!-- AYUDA: Presenten resultados vinculados con preguntas y objetivos. Incluyan
tablas o figuras con títulos, leyendas y archivos de origen. Describan aquí lo
obtenido; interprétenlo en Discusión.
EJEMPLO: Tabla de presencia y ausencia generada por src/compare_genes.py.
PRIMERA SESIÓN: dejen esta sección como pendiente. -->


1. **Tabla de expresión diferencial anotada.** De los 5394 genes evaluados, 4967 obtuvieron un valor de significancia corregido; los 427 restantes fueron apartados por el filtrado independiente del modelo por carecer de poder estadístico. Aplicando los umbrales de FDR < 0.05 y |log2FoldChange| > 1 se identificaron **1458 genes diferencialmente expresados: 830 inducidos y 628 reprimidos** durante la infección. Los resultados completos se encuentran en `results/tables/deg_results_completa.csv` y el subconjunto significativo en `deg_results_significativos.csv`. Para su interpretación visual se generaron un volcano plot y un mapa de calor:

	<img width="714" height="564" alt="imagen" src="https://github.com/user-attachments/assets/29a603cf-7f18-49e5-bb50-8e996eea29b0" />
	
	Figura 1. Volcano plot de los genes inducidos y reprimidos a partir del análisis de expresión diferencial

	<img width="1129" height="441" alt="imagen" src="https://github.com/user-attachments/assets/142c00d0-6d39-4e1a-90be-d90469c92eb6" />
	Figura 2. Heatmaps, de lado izquierdo se muestran los 20 genes más inducidos, de lado derecho los 20 genes más reprimidos

3. **Validación con controles biológicos.** Los cuatro controles positivos resultaron inducidos con significancia alta y en el orden de magnitud descrito en la literatura para la activación del regulón *hrp* durante la colonización del apoplasto:

| Locus | Gen | Función | log2FC | FDR | Clasificación |
|:--|:--|:--|--:|--:|:--|
| PSPTO2133 | *pvdS* | Adquisición de hierro | +6.86 | 1.2 × 10⁻⁶⁸ | Inducido |
| PSPTO4001 | *avrPto1* | Efector tipo III | +6.05 | 4.2 × 10⁻⁹⁶ | Inducido |
| PSPTO4709 | *cmaA* | Biosíntesis de coronatina | +5.68 | 4.9 × 10⁻¹⁸⁰ | Inducido |
| PSPTO1404 | *hrpL* | Regulador maestro del T3SS | +4.90 | 2.3 × 10⁻⁶¹ | Inducido |
| PSPTO1287 | *gap* | Gliceraldehído-3-fosfato deshidrogenasa | +1.44 | 6.6 × 10⁻⁹ | Inducido |
| PSPTO1745 | *gyrA* | ADN girasa, subunidad A | +0.99 | 1.8 × 10⁻⁷ | Sin cambio |

De los controles negativos, *gyrA* se mantuvo por debajo del umbral de cambio. *gap* lo superó ligeramente, con +1.44; dado que la gliceraldehído-3-fosfato deshidrogenasa participa en la glucólisis y su expresión depende del estado metabólico, no constituye un gen constitutivo en sentido estricto, por lo que un cambio moderado entre dos ambientes tan distintos es esperable.
 
3. **Enriquecimiento funcional de rutas KEGG.** A partir de la tabla de expresión diferencial se evaluó qué rutas biológicas están sobrerrepresentadas entre los genes que cambian su expresión, mediante un análisis de sobrerrepresentación (ORA) con prueba hipergeométrica y corrección por comparaciones múltiples. Entre los genes inducidos destacan rutas asociadas con la síntesis proteica (ribosoma), el metabolismo del carbono, la glucólisis/gluconeogénesis, el ciclo del ácido cítrico y la biosíntesis de metabolitos secundarios. Adicionalmente se identificaron 54 genes pertenecientes al Sistema de Secreción Tipo III (SST3), el principal mecanismo de virulencia de la bacteria. Los resultados se encuentran en `results/tables/ora_significant_pathways.csv` y la clasificación del SST3 en `degs_t3ss_clasificados.csv` y `resumen_t3ss.csv`:

   <img width="1396" height="615" alt="imagen" src="https://github.com/user-attachments/assets/82d2f5b7-e1f6-4c56-b3ef-c78392e2704c" />
   
   Figura 3. Rutas KEGG enriquecidas entre los genes inducidos durante la infección.
   
   <img width="1132" height="514" alt="imagen" src="https://github.com/user-attachments/assets/47d15963-db72-41b5-a466-70ddd0c736ca" />
   
   Figura 4. Rutas KEGG enriquecidas entre los genes reprimidos durante la infección.

4. **Genes sin función asignada como candidatos a virulencia.** De los 5394 genes evaluados, 525 están anotados únicamente como proteína hipotética, **122 presentan expresión diferencial** (31 inducidos y 91 reprimidos). Para priorizar candidatos se examinó el contexto genómico de los de mayor cambio, bajo el criterio de que en bacterias los genes contiguos suelen compartir regulación. El candidato principal es **PSPTO1385** (log2FoldChange = +6.32), la proteína hipotética más inducida de toda la tabla, que se localiza dentro de la isla de patogenicidad *hrp* flanqueada por genes estructurales del SST3 y se coinduce con ellos a magnitud equivalente. Un segundo candidato, PSPTO1372 (+5.40), se ubica entre un efector HopN1 y una chaperona CesT. Los resultados se encuentran en `results/tables/degs_funcion_desconocida.csv`:

| Locus | Producto anotado | log2FC |
|:--|:--|--:|
| PSPTO1383 | Subunidad de la varilla interna del SST3 | +6.56 |
| PSPTO1384 | Anillo de membrana interna del SST3 | +4.37 |
| **PSPTO1385** | **Proteína hipotética** | **+6.32** |
| PSPTO1386 | Proteína estatórica SctL del SST3 | +6.20 |
| PSPTO1387 | Proteína HrpF de secreción tipo III | +5.83 |

   Tabla 1. Contexto genómico del candidato PSPTO1385 y comportamiento de sus genes vecinos.

   <img width="1093" height="376" alt="imagen" src="https://github.com/user-attachments/assets/78c3d18c-d5bb-4698-90ee-dedead717cd5" />
   
   Figura 5. Ubicación genómica de PSPTO1385 dentro de la isla de patogenicidad *hrp*.


## 12. Discusión

<!-- AYUDA: Interpreten los resultados, expliquen si responden las preguntas,
compárenlos con los antecedentes y señalen limitaciones. No repitan únicamente
los valores.
EJEMPLO: La distribución observada sugiere..., aunque la interpretación está
limitada por la calidad de las anotaciones. -->
Los resultados del análisis transcriptómico diferencial de _Pseudomonas syringae pv. tomato DC3000_ permitieron explorar los cambios en la expresión génica asociados con la infección de _Arabidopsis thaliana_, así como identificar procesos biológicos y genes candidatos que podrían participar en la adaptación bacteriana al ambiente del hospedero.

***Cambios en la expresión génica durante la infección***.

En relación con nuestra primera pregunta de investigación, el análisis de expresión diferencial mediante InMoose permitió identificar genes cuya expresión varía significativamente entre las condiciones in planta e in vitro. La existencia de genes inducidos y reprimidos sugiere que la bacteria modifica su programa transcripcional en respuesta a las condiciones del ambiente vegetal.

Esta observación es consistente con lo descrito por Nobori et al. (2018), quienes demostraron que el transcriptoma de _P. syringae DC3000_ experimenta cambios durante su interacción con _A. thaliana_, evidenciando la importancia del ambiente del hospedero en la regulación de la expresión bacteriana [5]. Sin embargo, la separación observada entre condiciones mediante el análisis de componentes principales (PCA) no permite atribuir todas las diferencias exclusivamente al proceso de infección, ya que también pueden intervenir factores relacionados con las condiciones de cultivo, la variabilidad biológica y las características técnicas de las muestras.

***Funciones biológicas asociadas con la adaptación al hospedero***

Respecto a la segunda pregunta de investigación, el análisis de enriquecimiento funcional mediante KEGG permitió identificar procesos biológicos sobrerrepresentados entre los genes inducidos durante la infección, principalmente relacionados con la síntesis proteica, el metabolismo del carbono y la biosíntesis de metabolitos secundarios.
El enriquecimiento de la ruta del ribosoma podría reflejar cambios en la demanda de síntesis proteica durante la interacción bacteria-hospedero. Por otra parte, la representación significativa de rutas como glucólisis/gluconeogénesis, metabolismo del carbono y ciclo del ácido cítrico sugiere modificaciones en la expresión de genes asociados con la utilización de nutrientes y el metabolismo energético.

Estas observaciones son compatibles con los antecedentes descritos por Xin et al. (2018), quienes destacan que la capacidad de _P. syringae_ para colonizar el tejido vegetal depende, entre otros factores, de su adaptación a las condiciones fisiológicas y nutricionales del apoplasto [3]. Asimismo, el enriquecimiento de la biosíntesis de metabolitos secundarios podría relacionarse con cambios en procesos metabólicos relevantes para la fisiología bacteriana. No obstante, se requiere examinar los genes específicos asociados con esta ruta antes de proponer una relación directa con mecanismos de virulencia.

Es importante considerar que el enriquecimiento estadístico de una ruta no demuestra necesariamente su activación funcional, ya que el análisis se basa en la representación de genes diferencialmente expresados y no en mediciones directas de actividad enzimática, abundancia proteica o flujo metabólico.

***Genes con función desconocida y su posible relación con la patogenicidad***

La tercera pregunta de investigación se abordó mediante la identificación de genes diferencialmente expresados cuya función se encuentra anotada como hipotética o poco caracterizada. Entre los candidatos identificados destacan genes como PSPTO2142, PSPTO1210, PSPTO1385 y PSPTO1372, que presentaron una mayor expresión durante la infección. Estos resultados permiten proponerlos como candidatos para estudios posteriores, aunque su patrón de expresión no constituye evidencia suficiente para atribuirles una función específica en la patogenicidad.

Adicionalmente, el análisis del contexto genómico permitió examinar regiones que contienen genes relacionados con el sistema de secreción tipo III (T3SS), incluyendo efectores, proteínas estructurales y reguladores transcripcionales. Este sistema constituye uno de los principales mecanismos de virulencia de _P. syringae_, debido a su participación en la transferencia de proteínas efectoras hacia las células vegetales, como describen Buell et al. (2003) y Xin et al. (2018) [3,4].

La proximidad genómica de genes hipotéticos a regiones asociadas con la virulencia puede proporcionar información complementaria para priorizar candidatos. Sin embargo, esta asociación debe interpretarse con precaución, puesto que la cercanía entre genes no demuestra que compartan regulación, función o participación directa en el proceso de infección.

***Limitaciones e implicaciones del análisis***

Los resultados dependen de la calidad de los datos públicos y de las anotaciones funcionales disponibles. Aunque la integración con GenBank y KEGG permitió ampliar la interpretación biológica, algunos genes carecen de funciones caracterizadas o de correspondencia con rutas conocidas, lo que limita su clasificación funcional. En conjunto, los resultados responden a las preguntas de investigación desde una perspectiva bioinformática: permiten identificar genes diferencialmente expresados, reconocer procesos biológicos enriquecidos y seleccionar genes poco caracterizados para investigaciones posteriores.


## 13. Conclusiones

Con los resultados obtenidos del análisis de expresión diferencial y enriquecimiento funcional se cumplieron con los objetivos y preguntas de investigación.

A partir del análisis de expresión diferencial realizado con la herramienta Inmooose se encontraron 1458 genes, 830 inducidos y 628 reprimidos que responde a nuestra
primera pregunta de investigación.

Para conocer en que procesos estaba implicados, se realizó un enriquecimiento con KEGG  para asociar los genes cuyo patrón sea diferenciado del control con rutas biológicas, de los cuales solo 525 tenían anotación en la base de datos. Además, se encontraron 54 genes forman parte del Sistema de Secreción Tipo 3 (SST3), asociado a la patogenicidad, 
lo que contesta nuestra segunda pregunta de investigación, donde en la infección a la planta existe una mayor inducción de genes asociados a 
la patogenicidad. 

Se obtuvieron genes candidatos que podrían estar asociados a la virulencia por su patrón de expresión. Especialmente al analizar la ubicación cromosómica de PSPTO1385, se encontró en una región cercana con el SST3 que representa un logfoldchange de 6.32, el gen sin anotación más inducido. Lo que responde a nuestra 
tercera pregunta de investigación, sin embargo, es necesario realizar otros análisis en trabajos futuros para poder caracterizar de mejor manera su función 
y comprobar que es un gen asociado a la patogénesis que por ahora únicamente se cataloga como candidato, por ejemplo se podría hacer un BLASTN para 
conocer si tiene algún homólogo en otra cepa o bacteria y se podría elucidar su función o realizar una búsqueda más exhaustiva en otras bases de datos.


## 14. Disponibilidad, licencia y citación

<!-- AYUDA: Indiquen dónde está el código, bajo qué licencia puede reutilizarse
y cómo citarlo. Relacionen esta sección con LICENSE, CITATION.cff, codemeta.json,
release final y, cuando corresponda, un DOI.
EJEMPLO: Código en GitHub bajo MIT; cita disponible en CITATION.cff. -->

**Código:** [Código](https://github.com/ItayetziC/Plant_pathogen_interaction.git)
**Datos:** [Instrucciones de obtención](https://github.com/ItayetziC/Plant_pathogen_interaction.git)  
**Licencia del código:** [Licencia](https://github.com/ItayetziC/Plant_pathogen_interaction/blob/main/LICENSE)  
**Cómo citar:** [Para citar](https://github.com/ItayetziC/Plant_pathogen_interaction/blob/main/CITATION.cff)  
**Versión o release:** [Versión](https://github.com/ItayetziC/Plant_pathogen_interaction/releases/tag/v1.2.1)

## 15. Referencias

<!-- AYUDA: Registren publicaciones, datos, software y documentos consultados en
un formato uniforme. Incluyan DOI, URL o identificadores persistentes. Toda cita
del texto debe aparecer aquí.
EJEMPLO: Blattner, F. R. et al. (1997). The complete genome sequence of
Escherichia coli K-12. Science, 277(5331), 1453–1462.
https://doi.org/10.1126/science.277.5331.1453 -->

1. Shao, X., Tan, M., Xie, Y., Yao, C., Wang, T., Huang, H., Zhang, Y., Ding, Y., Liu, J., Han, L., Hua, C., Wang, X., & Deng, X. (2021). Integrated regulatory network in Pseudomonas syringae reveals dynamics of virulence. Cell Reports, 34(13), 108920. 10.1016/j.celrep.2021.108920
2. Helmann, T. C., Deutschbauer, A. M., & Lindow, S. E. (2019). Genome-wide identification of Pseudomonas syringae genes required for fitness during colonization of the leaf surface and apoplast. Proceedings of the National Academy of Sciences (PNAS), 116(38), 18900–18910. 10.1073/pnas.1908858116
3. Xin, X.-F., Kvitko, B., & He, S. Y. (2018). Pseudomonas syringae: what it takes to be a pathogen. Nature Reviews Microbiology, 16(5), 316–328.10.1038/nrmicro.2018.17
4. Buell, C. R., Joardar, V., Lindeberg, M., Selengut, J., Paulsen, I. T., Gwinn, M. L., Dodson, R. J., Deboy, R. T., Durkin, A. S., Kolonay, J. F., Madupu, R., Daugherty, S., Brinkac, L., Beanan, M. J., Haft, D. H., Nelson, W. C., Zafar, N., Zhou, L., Liu, J., Yuan, Q., Khouri, H., Fedorova, N., Tran, B., Russell, D., Berry, K., Utterback, T., Van Aken, S. E., Feldblyum, T. V., D'Ascenzo, M., Deng, W.-L., Ramos, A. R., Alfano, J. R., Chatterjee, A. K., Delaney, T. P., Lazarowitz, S. G., Martin, G. B., Schneider, D. J., Tang, X., Bender, C. L., White, O., Fraser, C. M., & Collmer, A. (2003). The complete genome sequence of the Arabidopsis and tomato pathogen Pseudomonas syringae pv. tomato DC3000. Proceedings of the National Academy of Sciences (PNAS), 100(17), 10181–10186. 10.1073/pnas.1733273100
5. Nobori T., Velasquez A. C., Wu J., Kvitko B. H., Kremer J. M., Wang Y., He S. Y. and Tsuda K.(2018). Transcriptome landscape of a bacterial pathogen under plant immunity. Proc Natl Acad Sci U S A 115(13): E3055-E3064. 10.1073/pnas.180052911
6. Colange, M., Appé, G., Meunier, L. _et al._ Differential expression analysis with inmoose, the integrated multi-omic open-source environment in Python. _BMC Bioinformatics_  26, 160 (2025). https://doi.org/10.1186/s12859-025-06180-7
7. Cock, P. A., et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422–1423. 10.1093/bioinformatics/btp163 
8. Waskom, M. L. (2021). seaborn: statistical data visualization. Journal of Open Source Software, 6(60), 3021. 10.21105/joss.03021
---

<!-- ORIENTACIÓN PARA LAS DOS PRIMERAS SESIONES:
Completen Información general, Resumen provisional, secciones 1 a 7, fuentes de
datos preliminares y plan de trabajo. Metodología, Resultados, Discusión,
Conclusiones Y Disponibilidad evolucionarán durante el
semestre. Sustituyan las indicaciones entre corchetes por contenido del equipo. -->
