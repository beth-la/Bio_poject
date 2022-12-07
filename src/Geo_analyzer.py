'''
Name
    Geo_analyzer.py
    
Version
    1.0
    
Authors
    Delgado Gutierrez Diana 
    Lopez Angeles Brenda E.
    Plasencia Perez Victor Ulises
    Ávila Silva Rogelio Lael
    
Descripcion
    
Category
    Biopython
    
Usage
    Python Geo_analyzer.py [-h] -o ORGANISM -f FEATURE
    
Arguments
    -h --help

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


# Obteniendo el query con la funcion entrez_query.
query = entrez_query(arguments.ORGANISM, arguments.FEATURE)
# Dependiendo del modo se ejecuta un código u otro.

if arguments.MODE == '1':
    if not arguments.ORGANISM or not arguments.FEATURE:
        print('\nFaltan argumentos, consultar opcion -h\n')
        exit(0)
    # Obteniendo los Ids de GSE asociados al query.
    gse_ids = format_gse(query)

    # Imprimir los IDs de GSE asociados
    print(" ".join(gse_ids))

elif arguments.MODE == '2':
    if not arguments.GEOid:
        print('\nEs necesario especificar un ID.\n')
        exit(0)
    gse_object_extract(arguments.GEOid)

else:
    print(
        '\nEs anecesario especificar un modo valido.\n')
