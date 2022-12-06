'''
MODULE NAME:
    Entrez module

DESCRIPTION:
    Modulo que contiene funciones que nos permiten trabajar con Entrez.

VERSION:
    1.0
    
AUTHORS
    Delgado Gutierrez Diana 
    Lopez Angeles Brenda E.
    Plasencia Perez Victor Ulises
    
Functions
    entrez_query
    format_gse
    
'''
# Librerias utilizadas por el modulo
from Bio import Entrez
import re 

def entrez_query(organism, feature):
    '''
    Genera un query que pueda ser utilizado con Entrez.
    Parameters:
        organism (str): Organismo para el cual se desea generar el query de busqueda.
        feature (str): Caracteristica que acompaña al organismo buscado 
    Returns:
        query (str): Query generado. 
    '''
    query= f"({feature}) AND ({organism} [Organism]) AND (Expression profiling by array [DataSet Type])" 
    return(query)

def format_gse(query):
    '''
    Obtiene una lista con los IDs GSE de la base de datos GEO.
    Parameters:
        query (str): Query para el cual se desea realizar la busqueda.
    Returns:
        gse_list (list): Lista de IDs listos para ser procesados por GEOparse. 
    '''
    Entrez.email = "blopez@lcg.unam.mx"
    handle = Entrez.esearch(db='gds', term=query, retmax=5, usehistory=True)
    result = Entrez.read(handle)
    handle.close()
    uid_regex = re.compile('[1-9]+0+([1-9]+[0-9]*)')
    gse_list = ['GSE' + uid_regex.match(uid).group(1) for uid in result['IdList']]
    return(gse_list)