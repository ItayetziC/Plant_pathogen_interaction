# Análisis transcriptómico diferencial de Pseudomonas syringae durante la infección de Arabidopsis thaliana

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
| Última actualización | [30/08/2026] |
| Repositorio | https://github.com/ItayetziC/Plant_pathogen_interaction |

<!-- AYUDA: Mantengan actualizados el estado, la fecha y la versión. La versión
debe coincidir con una etiqueta de Git cuando exista una entrega identificable. -->

## Resumen del proyecto

<!-- AYUDA: Expliquen en un párrafo el problema, su relevancia, los datos, la
solución propuesta y el resultado principal. Conviene actualizarlo al final.
EJEMPLO: Se analizarán genomas de E. coli disponibles en NCBI para identificar
y comparar genes de resistencia. Se desarrollará un flujo reproducible en Python
que obtenga los datos y genere tablas y visualizaciones comparativas. -->

[Redacten aquí el resumen.]

## 1. Contexto y antecedentes

<!-- AYUDA: Presenten la información necesaria para comprender el proyecto.
Definan conceptos biológicos y computacionales, describan qué se conoce y citen
trabajos, datos o herramientas relacionados.
PREGUNTAS: ¿Cuál es el fenómeno de interés? ¿Qué debe conocer quien lea el
reporte? ¿Qué métodos o herramientas se han utilizado anteriormente?
EJEMPLO: La resistencia antimicrobiana es un problema de salud pública. Aunque
existen bases especializadas, comparar varios genomas requiere integrar datos
procedentes de distintos archivos. -->

_Pseudomonas syringae_ pv.tomato DC3000 es una bacteria Gram-negativa patógena que infecta 
tanto al tomate como a la planta modelo Arabidopsis thaliana, lo que la ha convertido en
un sistema de referencia más estudiados para descifrar los mecanismos moleculares de la 
patogénesis bacteriana y la inmunidad vegetal [1,2].

La bacteria para infectar a la planta ingresa a la hoja por los estomas, suprime la 
inmunidad vegetal para poder multiplicarse en el apoplasto y convertirlo en acuoso 
(inundación del espacio intercelular) para absorber nutrientes [2,3]. 

El control transcripcional del proceso de infección ha sido ampliamente estudiado para 
comprender la red de genes virulentos, descritos gracias a la secuenciación completa del 
genoma de DC3000 que también reveló un cromosoma de 6.4 Mb y dos plásmidos que codifican 
5,763 marcos de lectura abierta, entre ellos 298 genes de virulencia establecidos y 
putativos [4]. 

El programa transcripcional que la bacteria despliega realmente dentro del hospedero permanece
menos caracterizado ya que en tejido vegetal infectado el ARN bacteriano representa una fracción 
mínima del ARN total, dominado por los transcritos de la planta, sin embargo, en 2018 Nobori _et al,._ [5] desarrollaron métodos de enriquecimiento que permitieron superar este obstáculo y 
generaron el transcriptoma de DC3000 durante la infección de _A. thaliana_, depositado públicamente como `GSE103442`.

En este proyecto se realizará un análisis de expresión génica diferencial utilizando el 
conjunto de datasets de Buell _et al.,_ (2003) y Nobori _et al.,_ (2018) para comparar 
estadísticamente la cantidad de ARN producido por los genes, es decir su expresión, en dos 
condiciones distintas: durante la infección a _A. thaliana_ y en condiciones no infectivas 
(crecimiento en medio de cultivo) [4,5]. 

Apriximadamente 658 estaban anotados como proteínas hipotéticas, lo que sugiere que una cantidad 
considerable de genes durante la infección carece de función asignada [5], de modo que se podrían identificar cambios de expresión cuyo significado biológico permanece indeterminado.

El análisis se implementará en Python con PyDESeq2 [6], y la anotación funcional de los genes 
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
importancia económica mundial para la cual aún no existen medidas de control efectivas (4). 
Diseñarlas requiere saber qué genes utiliza la bacteria durante la infección real, sin embargo, hasta ahora las investigaciones sobre regulación de su virulencia proviene de experimentos en medios de cultivo, es decir, no reproducen las condiciones del apoplasto. Además solo el 61% de sus 5763 genes tiene función asignada, por lo que son necesarios más anáisis para comprender su mecanismo completo de patogenicidad.


## 3. Justificación

<!-- AYUDA: Expliquen por qué vale la pena realizar el proyecto, cuál es su
relevancia biológica, científica, técnica o social y quién podría beneficiarse.
EJEMPLO: Un flujo automatizado reducirá errores y permitirá repetir el análisis
con los mismos datos, parámetros y versiones del software. -->

_Pseudomonas_ es la causante de más de 50 tipos de enfermedades en diferentes cultivos que afectan a nivel global tanto económicamente como socialmente. 
Aún no existe una estrategia biotecnológica para controlar a _Pseudomonas syringae_ pv. tomato DC3000 ni se conocen todos los genes involucrados en su patogenicidad que podrían ser clave para su desarrollo. 
Aunque los datos transcriptómicos necesarios para reducir este vacío ya son públicos, aún permanecen subutilizados y determinar cuáles cambian de forma significativa podría elucidar genes candidatos que se prioricen para estudios funcionales que cambien dependiendo de su expresión al infectar a la planta modelo _A. thaliana_.

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

**Pregunta:** ¿Qué genes de _Pseudomonas syringae_ pv.tomato DC3000 presentan cambios 
significativos en su expresión durante la infección de _Arabidopsis thaliana_ en comparación 
con condiciones de crecimiento no infectadas?  

**Evidencia necesaria:** Matriz de expresión génica y réplicas experimentales de ambas condiciones,
a partir de las cuales se obtenga una tabla de genes diferencialmente expresados con sus valores 
y significancia. Los resultados podrán visualizarse mediante un volcano plot y un heatmap.


### Pregunta 2

**Pregunta:** ¿Qué funciones biológicas y procesos asociados con la adaptación al hospedero y la 
patogénesis están representados entre los genes diferencialmente expresados durante la infección?

**Evidencia necesaria:** Anotación funcional de los genes diferencialmente expresados y clasificación de estos según sus funciones, particularmente aquellas relacionadas con secreción y efectores, transporte y adquisición de nutrientes, regulación, adhesión, respuesta al estrés y otros mecanismos asociados con la virulencia.

### Pregunta 3

**Pregunta:** ¿Existen genes diferencialmente expresados durante la infección cuya función se 
encuentre poco caracterizada o sea desconocida y que, por su patrón de expresión, puedan 
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

- **Análisis bioinformático** (Control de calidad de datos, cuantificación de la expresión génica)
- **Muestras y datos:** Datos secundarios de RNA-seq (Serie GEO `GSE103442` / BioProject `PRJNA419365`, 
Nobori et al., 2018) que abarcan condiciones de infección (_in planta_) frente a controles de 
crecimiento no infectivo (_in vitro_).
- **Genoma de referencia:** Ensamblado de RefSeq `GCF_000007805.1` para el alineamiento y anotación 
funcional de la bacteria.

### No incluye

- **Validación experimental** (Ensayos de laboratorio (_wet-lab_), como RT-qPCR, mutagénesis ni infecciones 
_in vivo_).
- **Otras cepas u organismos:** No se evaluarán otras cepas de _P. syringae_ ni otros patógenos u 
hospederos alternativos.    
-   **Transcriptómica del hospedero:** No se caracterizará la respuesta transcripcional de 
_Arabidopsis thaliana_, enfocándose únicamente en el transcriptoma bacteriano.
-   **Secuenciación de nuevo (_de novo_):** No se generarán datos de secuenciación propios ni se 
realizará ensamblado _de novo_.

### Limitaciones conocidas

- **Origen de los datos:** Dependencia total de la calidad, profundidad de secuenciación y diseño 
experimental de los datos públicos depositados en NCBI por Nobori _et al,._ (2018) y Buell _et al,._ (2003).
- **Falta de anotación certera:** A pesar de identificar cambios de expresión, no se podrá asegurar la función de genes putativos en el proceso de infección.

## 7. Propuesta de solución

<!-- AYUDA: Describan el producto o estrategia que podría resolver el problema.
Es una propuesta inicial y puede cambiar. Expliquen sus componentes, no sólo las
tecnologías.
EJEMPLO: Un programa modular recibirá identificadores, descargará archivos,
extraerá genes, almacenará resultados y generará visualizaciones. -->

[Describan aquí la solución propuesta.]

### 7.1 Resultado o producto esperado

<!-- AYUDA: Indiquen el entregable concreto: programa, paquete, flujo de análisis,
base de datos, visualizaciones u otro producto.
EJEMPLO: Repositorio ejecutable con scripts, datos de prueba, documentación,
tabla comparativa y figuras regenerables. -->

[Describan aquí el producto.]


## 8. Datos

### 8.1 Fuentes de datos

<!-- AYUDA: Incluyan institución, base de datos, URL, identificador, versión o
fecha de consulta y condiciones de uso. No todos los proyectos usarán NCBI.
EJEMPLO: NCBI RefSeq, GCF_000005845.2, consultado el dd/mm/aaaa. -->

| Fuente | Identificador o versión | URL | Fecha de consulta | Licencia o condiciones |
|:--|:--|:--|:--|:--|
| NCBI Assembly / RefSeq (NCBI / NIH) | Genoma de referencia _Pseudomonas syringae_ pv. _tomato_ DC3000	`GCF_000007805.1` (NC_004578.1) |	[https://www.ncbi.nlm.nih.gov/assembly/29388](https://www.ncbi.nlm.nih.gov/assembly/29388)	|	30/08/2026	|	Dominio Público. Uso libre para investigación académica y comercial.
| NCBI GEO (Gene Expression Omnibus)	|	Datos de RNA-seq de _P. syringae_ DC3000 (_in planta_ vs. _in vitro_) Serie GEO: `GSE103442` (BioProject: `PRJNA401390`)	|	https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103442	|	30/08/2026	|	Datos de acceso abierto. Reutilización libre citando la publicación de Nobori et al. (2018).

### 8.2 Características de los datos

<!-- AYUDA: Describan organismos, muestras, variables, formatos, versiones, tamaño y otros
atributos necesarios para interpretar los datos.
EJEMPLO: Archivos FASTA y GFF3 de 20 genomas completos de E. coli. -->

**1. Genoma de referencia de _Pseudomonas syringae_**
-   **Organismo:** _Pseudomonas syringae_ pv. _tomato_ (cepa DC3000).
-   **Muestras / Tipo de datos:** Genoma bacteriano completo de referencia y su anotación funcional. 
-   **Variables / Atributos:** Secuencia completa de nucleótidos del cromosoma bacteriano 
(~6.39 Mb, ~5,600 genes codificantes de proteínas) y sus plásmidos asociados (pDC3000A y pDC3000B). 
Incluye coordenadas genómicas de exones, CDS, rRNAs, tRNAs y nombres de genes.  
-   **Formatos:** Archivo de secuencia en formato **FASTA** (`.fna`) y archivo de anotación genómica 
en formato **GFF3** (`.gff`) o **GTF**.
-   **Versión:** Ensamblado de RefSeq `GCF_000007805.1` (Secuencia de acceso al cromosoma principal: 
`NC_004578.1`).
-   **Tamaño:** ~2.0 MB comprimido (`.gz`).

**2.  Descripción de los datos transcriptómicos de _Pseudomonas syringae_**

-   **Organismo:** _Pseudomonas syringae_ pv. _tomato_ (cepa DC3000).
-   **Muestras:** Muestras de ARN bacteriano extraído directamente de dos condiciones experimentales:
    1.  _In planta:_ Bacterias aisladas o recolectadas del apoplasto de _Arabidopsis thaliana_ 
    durante el proceso de infección.
    2.  _In vitro:_ Bacterias cultivadas en medio sintético de laboratorio (condición control/no infectiva).
-   **Variables / Atributos:** Niveles de abundancia de transcritos de ARN (expresión génica bacteriana) 
en respuesta al entorno del hospedero. Mide la cantidad de lecturas mapeadas por gen bacteriano para 
comparar la activación de factores de virulencia, efectores y metabolismo adaptativo.
-   **Formatos:**
    -   **Datos procesados:** Matriz de conteos (_read counts table_) en texto plano delimitado 
    por tabuladores (`.txt` o `.tsv` comprimido en `.gz`), con filas representando IDs de genes de 
    DC3000 y columnas representando las réplicas biológicas.
    -   **Datos brutos:** Archivos de secuenciación de alto rendimiento (_Illumina_) en formato 
    **FASTQ** comprimido (`.fastq.gz`) almacenados en SRA.
-   **Versión:** Serie GEO `GSE103442` (Estudio de Nobori et al., 2018).
-   **Tamaño:** ~5–15 MB para la tabla de datos procesados/conteos; varios Gigabytes (~10–30 GB) 
para los archivos brutos FASTQ completos.


### 8.3 Organización de los datos

<!-- AYUDA: Muestren la estructura prevista. No suban datos sensibles, tokens,
contraseñas ni archivos grandes. Usen .gitignore y documenten cómo obtener lo
que no se guarde en Git.
EJEMPLO: data/raw conserva originales y data/processed los derivados. -->

```text
proyecto/
├── data/
│   ├── raw/
│   └── processed/
├── docs/
├── notebooks/
├── results/
├── src/
└── tests/
```

### 8.4 Diccionario o formato de los datos

<!-- AYUDA: Describan campos o columnas relevantes. Incluyan fragmentos pequeños
cuando ayuden a comprender el formato, pero no archivos completos.
EJEMPLO: En GFF3, seqid identifica la secuencia; type indica gene, CDS, etc. -->

| Archivo o conjunto | Campo/columna | Tipo | Descripción | Valores o unidades |
|:--|:--|:--|:--|:--|
| [Archivo] | [Campo] | [Tipo] | [Descripción] | [Valores] |

## 9. Metodología

<!-- AYUDA: Esta sección evolucionará. Primero describan el plan y después
actualícenla con lo que realmente ejecutaron, incluidos parámetros y decisiones. -->

### 9.1 Etapas del análisis o desarrollo

<!-- AYUDA: Describan la secuencia desde la obtención de datos hasta la validación
de resultados. Relacionen cada etapa con una pregunta u objetivo.
EJEMPLO: descarga, validación, transformación, análisis, visualización y pruebas. -->

1. [Etapa 1]
2. [Etapa 2]
3. [Etapa 3]

### 9.2 Herramientas y tecnologías

<!-- AYUDA: Registren lenguajes, bibliotecas y programas con sus versiones y
propósito. No incluyan credenciales.
EJEMPLO: Python 3.x; Biopython para leer formatos biológicos; Seaborn para
visualización. -->

| Herramienta | Versión | Propósito |
|:--|:--|:--|
| [Herramienta] | [Versión] | [Uso] |


### 9.3 Estrategia de validación

<!-- AYUDA: Expliquen cómo comprobarán código y resultados: pruebas unitarias,
datos conocidos, comparación con otra herramienta o revisión manual.
EJEMPLO: Se compararán cinco anotaciones conocidas y se probarán entradas
válidas, identificadores inexistentes y archivos incompletos. -->

[Describan aquí la validación.]

## 10. Plan de trabajo


### 10.1 Distribución de responsabilidades

<!-- AYUDA: Definan responsabilidades iniciales sin aislar a cada integrante.
Toda contribución importante debe ser revisada mediante Pull Request por otra
persona.
EJEMPLO: Ana desarrolla la descarga y revisa el módulo de visualización. -->

| Integrante | Responsabilidad principal | Responsabilidad de revisión |
|:--|:--|:--|
| [Nombre] | [Responsabilidad] | [Qué o a quién revisará] |

### 10.2 Riesgos y alternativas

<!-- AYUDA: Identifiquen situaciones que podrían impedir o retrasar el proyecto
y definan una alternativa.
EJEMPLO: Los datos requieren demasiado almacenamiento; alternativa: reducir el
número de genomas usando criterios documentados. -->

| Riesgo | Probabilidad | Impacto | Prevención o alternativa |
|:--|:--|:--|:--|
| [Riesgo] | Baja/Media/Alta | Bajo/Medio/Alto | [Acción] |

## 11. Resultados

<!-- AYUDA: Presenten resultados vinculados con preguntas y objetivos. Incluyan
tablas o figuras con títulos, leyendas y archivos de origen. Describan aquí lo
obtenido; interprétenlo en Discusión.
EJEMPLO: Tabla de presencia y ausencia generada por src/compare_genes.py.
PRIMERA SESIÓN: dejen esta sección como pendiente. -->

> Estado: pendiente. Se completará durante el desarrollo.



## 12. Discusión

<!-- AYUDA: Interpreten los resultados, expliquen si responden las preguntas,
compárenlos con los antecedentes y señalen limitaciones. No repitan únicamente
los valores.
EJEMPLO: La distribución observada sugiere..., aunque la interpretación está
limitada por la calidad de las anotaciones. -->

> Estado: pendiente. Se completará después de obtener resultados.

## 13. Conclusiones

<!-- AYUDA: Sinteticen qué se aprendió, qué preguntas se respondieron y si se
alcanzaron los objetivos. Incluyan aportes, limitaciones y trabajo futuro. No
introduzcan resultados nuevos.
EJEMPLO: El flujo permitió identificar..., pero será necesario incorporar... -->

> Estado: pendiente. Se completará al finalizar el proyecto.


## 14. Disponibilidad, licencia y citación

<!-- AYUDA: Indiquen dónde está el código, bajo qué licencia puede reutilizarse
y cómo citarlo. Relacionen esta sección con LICENSE, CITATION.cff, codemeta.json,
release final y, cuando corresponda, un DOI.
EJEMPLO: Código en GitHub bajo MIT; cita disponible en CITATION.cff. -->

**Código:** [URL]  
**Datos:** [URL, identificador o instrucciones]  
**Licencia del código:** [Licencia]  
**Cómo citar:** [Referencia o enlace a CITATION.cff]  
**Versión o release:** [URL]

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
6. Muzellec, B., Telenczuk, M., Cabeli, V., & Andreux, M. (2023). PyDESeq2: a python package for bulk RNA-seq differential expression analysis. Bioinformatics, 39(9), btad547. 10.1093/bioinformatics/btad547
7. Cock, P. A., et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422–1423. 10.1093/bioinformatics/btp163 
8. Waskom, M. L. (2021). seaborn: statistical data visualization. Journal of Open Source Software, 6(60), 3021. 10.21105/joss.03021
---

<!-- ORIENTACIÓN PARA LAS DOS PRIMERAS SESIONES:
Completen Información general, Resumen provisional, secciones 1 a 7, fuentes de
datos preliminares y plan de trabajo. Metodología, Resultados, Discusión,
Conclusiones Y Disponibilidad evolucionarán durante el
semestre. Sustituyan las indicaciones entre corchetes por contenido del equipo. -->
