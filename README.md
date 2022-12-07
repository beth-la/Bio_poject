# Bio_poject

## Resultados
entrez_module.py <- Las funciones contenidas nos permiten trabajar con la búsqueda de datos de expresión a través de GEOparse.
entez_module nos produce el query generado por el usuario y una lista de IDs listos para GEOparse.

Geo_analyzer.py <- Imprime en pantalla título, resumen, tipo, IDsde PubMed y IDs de plataforma del GSE producido en entrez_module.py
Regresa un data frame con los genes diferencialmente expresados con base a los datos de la plaraforma de cada GSE. 

dexs_class.py <- Hace gráficas con los datos de expresión diferencial de tabla que regresa de Geo_analyzer.py. Tiene dos atributos: un cluster map de los datos de expresión diferencial de controles y los tratados y otro que los IDs de los genes.

## Conclusion
Se diseñó una herramienta que ayuda a buscar diseños experimentales con datos de expresión obtenida mediante microarreglos que permitan analizar a través de gráficas y tablas la expresión diferencias de los mismos datos. A través de una cadena el usuario puede buscar datos de expresión en GEO que le permiten explorar y hacer análisis sobre esos datos para su posterior intepretación.
Esta herramienta es un primer acercamiento para saber los resultados esperados antes de tener un apropuesta experimental. Asimismo nos permite obtener una base para el análisis de datos nuevos y cómo obtenerlos a través de una metodología reproducible.
Tambien es una inicialización hacia nuevos proyectos y nuevas preguntas que se espera sean resuletas con un análisis más profundo basándose en los datos proporcionados con este.