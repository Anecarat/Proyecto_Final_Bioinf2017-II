# README
## Proyecto final de Bioinformática: Análisis de la comunidad de bacterias asociadas a las raíces y hojas de *Boecheria stricta*

El presente proyecto contiene los datos, scripts e información necesaria para realizar, desde el pre-procesamiento de datos crudos, hasta los análisis de diversidad para secuencias del gen del 16S del rRNA obtenidas de comunidades de bacterias asociadas a *B. stricta*. 

Los datos fueron tomados del [repositorio](http://dx.doi.org/10.1038/ncomms12151) de datos de Wagner et al. (2016) y los scripts generados están basados en los scripts de este mismo repositorio y el [workflow](https://f1000research.com/articles/5-1492/v1) de bioconductor de Callaghan et al.

## Contenido

En esta carpeta encontrarás los siguientes archivos: 
- **0.MasterScript:** En este script se encuentra la información de la sesión de R con la que se trabajó el proyecto, además de los comandos necesarios para crear las carpetas, acomodar los archivos e instalar todas las paqueterías necesarias para correr los análisis. 
- **Initial:** En esta carpeta se encuentran los scripts y datos iniciales para realizar los análisis. 

## ¿Qué hace este proyecto?

Dentro de los archivos que se encuentran en la carpeta **Initial**, se encuentran otros cuatro scripts que realizan lo siguiente: 
- **1.PreProcessing:** Este script contiene las instrucciones necesarias para procesar las secuencias crudas, tal cual se entregan de la plataforma de IlluminaMiSeq. Incluye la eliminación de los barcodes y primers, análisis de calidad de las secuencias y el trimming de éstas, la unión de las muestras en una misma tabla, la eliminación de quimeras y generación de unidades taxonómicas operacionales (OTU).
- **2.AssignTaxonomy:** En este script se clasificarán taxonómicamente cada uno de los OTUs, se construirá el árbol filogenético de las muestras y se anexarán datos ambientales relevantes de las muestras a la tabla de OTUs.
- **3.DiversityAnalyses:** En este, se encuentran los comandos para calcular diferentes índices de diversidad $\alfa$ y $\beta$, además de curvas de rarefacción de cada una de las muestras. 
- **4.MultivariateAnalyses:** Finalmente, en este script se realizarán algunas pruebas multivariadas (aún por definir cuáles; no sé si son diferentes aproximaciones a lo mismo o cada una aporta informacón diferente) para analizar la similitud entre las muestras. 

