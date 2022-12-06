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