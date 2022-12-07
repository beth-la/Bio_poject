# Bio_poject

## Diseño

Se creó un programa con dos modos de uso principales.
En su primer modo de uso, realiza una búsqueda en la base de datos de Gene Expression Omnibus (GEO) con el módulo Entrez para la obtención de IDs de datos de expresión obtenidos mediante microarreglos de un organismo y una característica, ambas dadas por el usuario en los argumentos del programa.

Dentro del segundo modo de uso, un ID de GEO se toma para realizar un análisis de expresión diferencial con los datos de expresión asociados al ID, que regresa una tabla con los datos resultantes del análisis como los niveles de expresión, símbolo del gen, anotación y su ID de Ilumina. Con este fin, GEOparser es un módulo específico para el manejo de datos de microarreglos. Otros dos módulos, pandas y numpy, se emplearon para los cálculos necesarios para la determinación de genes con expresión diferencial y finalmente se diseñó una clase con el formato de la tabla de los genes resultantes con atributos que faciliten la visualización de los datos gráficamente.

## Autores

[Victor Plascencia](https://github.com/ulisesplaper), [Diana Delgado](https://github.com/dianadg159/python_class), [Brenda López](https://github.com/beth-la), [Rogelio Ávila](https://github.com/Roglavsil)

## Objetivo

Creación de un programa útil para análisis de expresión diferencial de datos de expresión de la base de datos Gene Expression Omnibus (GEO).

## Versionado

Las versiones se generaron manualmente

## Estructura

- data: Archivos input para programas

- results: Archivos output de programas

- src: Scripts de Python

- readme

## Resultados

## Conclusión

## Pre-requisitos

python 3.10.2
