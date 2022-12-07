'''
Name
    Geo_analyzer.py
    
Version
    1.2
    
Authors
    Delgado Gutierrez Diana 
    Lopez Angeles Brenda E.
    Plasencia Perez Victor Ulises
    Ávila Silva Rogelio Lael
    
Descripcion
    Modo 1: obtiene IDs de la db GEO asociados 
        a un estado(feature) y un organismo
    Modo 2: Regresa un DF con los DEGs y un heatmap
        dados un ID, una lista de IDs de controles, 
        una lista de IDs de muestras y un nombre
        de la columna de la tabla de la plataforma
        que contiene los datos de anotacion
Category
    Biopython
Usage
    Modo 1: Python Geo_analyzer.py -o ORGANISM -f FEATURE
    Modo 2: Python Geo_analyzer.py -ID GEOID -c CONTROLSID
        -s SAMPLESID -an ANNOTNAME
    
Arguments
     -h, --help            show this help message and exit
    -o ORGANISM, --ORGANISM ORGANISM
                        Organismo para el cual se desea encontrar información
    -f FEATURE, --FEATURE FEATURE
                        Caracteristica asociada al organismo
    -id GEOID, --GEOid GEOID
                        ID de la base de datos GEO
    -m MODE, --MODE MODE  Modo de uso del programa:

                        1: Obtener los ids de GEO asociados al término.
                        Argumentos requeridos: --ORGANISM, --FEATURE.

                        2: Análisis de expresión de un ID.
                        Argumentos requeridos: --GEOid
    -lfc LOGFOLDCHANGE, --logFoldChange LOGFOLDCHANGE
                        Logaritmo de duplicacion/reduccion de expresion
    -c CONTROLSID, --controlsID CONTROLSID
                        lista con los IDs de los controles a ser usados
    -s SAMPLESID, --samplesID SAMPLESID
                        lista con los IDs de las muestras a ser usadas
    -an ANNOTNAME, --annotName ANNOTNAME
                        nombre de la columna con los datos de anotacion

See also
    None
'''

# Importando modulos:
from entrez_module import entrez_query
from entrez_module import format_gse
import argparse
from Bio import Entrez
import GEOparse
from argparse import RawTextHelpFormatter
import numpy as np
import pandas as pd 


# Argumentos
arg_parser = argparse.ArgumentParser(
    description="", formatter_class=RawTextHelpFormatter)

arg_parser.add_argument("-o", "--ORGANISM",
                        help="Organismo para el cual se desea encontrar información",
                        required=False)

arg_parser.add_argument("-f", "--FEATURE",
                        help="Caracteristica asociada al organismo",
                        required=False)

arg_parser.add_argument("-id", "--GEOid",
                        help="ID de la base de datos GEO",
                        required=False)

arg_parser.add_argument("-m", "--MODE",
                        help="Modo de uso del programa:\n \n1: Obtener los ids de GEO asociados al término. \nArgumentos requeridos: --ORGANISM, --FEATURE.\n\n2: Análisis de expresión de un ID. \nArgumentos requeridos: --GEOid",
                        required=True)
arg_parser.add_argument("-lfc", "--logFoldChange",
                        help="Logaritmo de duplicacion/reduccion de expresion",
                        required=False, default=2)
arg_parser.add_argument("-c", "--controlsID",
                        help="lista con los IDs de los controles a ser usados",
                        required=False)
arg_parser.add_argument("-s", "--samplesID",
                        help="lista con los IDs de las muestras a ser usadas",
                        required=False)
arg_parser.add_argument("-an", "--annotName",
                        help="nombre de la columna con los datos de anotacion",
                        required=False)
arguments = arg_parser.parse_args()


def gse_object_extract(id):
    """
    Funcion: 
        Busta metadatos de experimentos con los GEO IDs proporcionados. 

    Args:
        ids (list): lista con GEO IDs
    Returns:
        object: class GEOparse.GEOTypes.GSE(name, metadata)
        A pantalla: Título, resumen, tipo y platform_id de los 
                    GEO IDs proporcionados.
    """
    # obtener el objeto gse
    gse = GEOparse.get_GEO(geo=id, destdir="./data/")

    # Title
    if 'title' in gse.metadata:
        print('Title: %s \n' % gse.metadata['title'][0])

    # Summary
    if 'summary' in gse.metadata:
        print('Summary: %s \n' % gse.metadata['summary'][0])

    # Type
    if 'type' in gse.metadata:
        print('Type: %s \n' % gse.metadata['type'][0])

    # PubMed IDs
    if 'pubmed_id' in gse.metadata:
        print('PubMed ID: %s \n' % gse.metadata['pubmed_id'][0])

    # platform_ids
    if 'platform_id' in gse.metadata:
        print('platform_id: %s \n' % gse.metadata['platform_id'][0])
    return (gse)

def make_DifExp_analysis(gse_objet, lfc, controls, samples, annot_column_name):
    """
    Funcion: 
        Obtiene los DEGs y su anotacion a partir de un gse_objet 

    Args:
        ids (list): lista con GEO IDs
        lfc (float): logaritmo del cambio de expresion
            en base 2.
        controls (list): lista con los IDs de los controles
        samples (list): lista con los IDs de las muestras
        annot_column_name(str): cadena con el nombre de la
            columna que contiene los nombres de los genes en 
            la plataforma
    Returns:
        DF con los DEGs, sus IDs, su anotacion y sus niveles de expr.
    """
    # Casting del lfc
    lfc = float(lfc)
    # Obtener la plataforma
    for platform in gse.gpls:
        break
    # Obtener los datos de expresion
    samples.extend(controls)
    pivoted_samples = gse.pivot_samples('VALUE')[samples]
    # Obtener el promedio por filas
    pivoted_samples_average = pivoted_samples.mean(axis=1)
    
    # Eliminamos el 25 % de datos mas bajos
    expression_threshold = pivoted_samples_average.quantile(0.25)
    # Obtenemos los ids asociadas a las muestras de expresion que 
    # se encuentran por arriba del valor del cuantil 
    expressed_probes = pivoted_samples_average[pivoted_samples_average >=expression_threshold].index.tolist()
    # Acceder a las muestras con los indices recuperados
    samples = pivoted_samples.loc[expressed_probes]
    
    # Obtener el disenio experimental
    experiments = {}
    for i,idx in enumerate(samples):
        tmp = {}
        tmp["Experiment"] = idx
        tmp["Type"] = "control" if idx in controls else "samples"
        experiments[i] = tmp
    experiments = pd.DataFrame(experiments).T
    
    # Se crea un dataframe con agrupando por muestras y controles.
    lfc_results = {}
    for tipo, group in experiments.groupby("Type"):
        lfc_results[tipo] = (samples.loc[:,group.Experiment].mean(axis = 1))
    lfc_results = pd.DataFrame(lfc_results)
    
    # Nombre de la columna que contiene los datos de anot
    interest_column = annot_column_name
    # Anotar con GLP
    lfc_result_annotated = lfc_results.reset_index().merge(gse.gpls[platform].table[["ID",interest_column]],left_on='ID_REF', right_on="ID").set_index('ID_REF')
    
    del lfc_result_annotated["ID"]
    
    # Remover celdas sin Entrez
    lfc_result_annotated = lfc_result_annotated.dropna(subset=[interest_column])
    # Remueve celdas con mas de un gene asociado
    lfc_result_annotated = lfc_result_annotated[~lfc_result_annotated.loc[:,interest_column].str.contains("///")]
    # Promediar los genes por cada celda
    lfc_result_annotated = lfc_result_annotated.groupby(interest_column).mean()
    # Obtner el logaritmo de base
    lfc_result_annotated = np.log2(lfc_result_annotated)

    # Obtener los genes diferencialmente expresados.
    DEGs = abs(lfc_result_annotated.control -lfc_result_annotated.samples) > lfc

    lfc_result_annotated.loc[:,'Diferentes'] = DEGs
    lcf_relevant = lfc_result_annotated.loc[lfc_result_annotated.Diferentes]
    # Aqui podriamos aniadir una excepcion
    lcf_relevant = lcf_relevant.reset_index().merge(gse.gpls[platform].table[[interest_column,'ID']], on=interest_column)
    lcf_relevant = lcf_relevant.set_index('ID')
    return(lcf_relevant)


# Obteniendo el query con la funcion entrez_query.
query = entrez_query(arguments.ORGANISM, arguments.FEATURE)
# Dependiendo del modo se ejecuta un código u otro.

if arguments.MODE == '1':
    # Verificar que se proporcionaron los 
    # argumentos necesarios
    if not arguments.ORGANISM or not arguments.FEATURE:
        print('\nFaltan argumentos, consultar opcion -h\n')
        exit(0)
    # Obteniendo los Ids de GSE asociados al query.
    gse_ids = format_gse(query)

    # Imprimir los IDs de GSE asociados
    print(" ".join(gse_ids))

elif arguments.MODE == '2':
    # Verificar que se proporcionaron los 
    # argumentos necesarios
    if not arguments.GEOid or not arguments.controlsID\
         or not arguments.samplesID or not arguments.annotName:
        print('\nFaltan argumentos para la opcion especificada.\n')
        exit(0)
    # Convertir los parametros a los TD
    # requeridos para la funcion
    gse = gse_object_extract(arguments.GEOid)
    samples = arguments.samplesID
    controls = arguments.controlsID
    samples_list = samples.split(',')
    controls_list = controls.split(',')
    # Llamar a la funcion para hacer 
    # el analisis de expresion diff
    print(make_DifExp_analysis(gse_objet= gse, lfc= arguments.logFoldChange, controls= controls_list, samples=samples_list, annot_column_name= arguments.annotName))

else:
    print(
        '\nEs anecesario especificar un modo valido.\n')
