'''
Name
    Geo_analyzer.py
    
Version
    1.0
    
Authors
    Delgado Gutierrez Diana 
    Lopez Angeles Brenda E.
    Plasencia Perez Victor Ulises
    
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

# Argumentos
arg_parser = argparse.ArgumentParser(description="")

arg_parser.add_argument("-o", "--ORGANISM",
                    help="Organismo para el cual se desea encontrar informacaion",
                    required=True)
                    
arg_parser.add_argument("-f", "--FEATURE",
                    help="Caracteristica asociada al organismo",
                    required=True)   

arg_parser.add_argument("-p","--PRINT",
                    help = "Imprimir a pantalla los IDs de GSE",
                    action = 'store_true',
                    required= False)         
arguments = arg_parser.parse_args()

# Obteniendo el query con la funcion entrez_query.
query= entrez_query(arguments.ORGANISM, arguments.FEATURE)

# Obteniendo los Ids de GSE asociados al query. 
gse_ids = format_gse(query)

# Si el usuario quiere imprimir los ids a pantalla: 
if arguments.PRINT:
    print(" ".join(gse_ids))
    

def gse_object_extract(ids):
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
    ids = ids[:2]
    gses_list = []
    for id in ids:
        gse = GEOparse.get_GEO(geo=id, destdir="./data/")
        # get accession ID
        acc = gse.get_accession()
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

        gses_list.append(gse)
    return(gses_list)