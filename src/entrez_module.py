from Bio import Entrez
import re 

def entrez_query(organism, feature):
    query= f"({feature}) AND ({organism} [Organism]) AND (Expression profiling by array [DataSet Type])" 
    return(query)

def format_gse(query):
    Entrez.email = "blopez@lcg.unam.mx"
    handle = Entrez.esearch(db='gds', term=query, retmax=5, usehistory=True)
    result = Entrez.read(handle)
    handle.close()
    uid_regex = re.compile('[1-9]+0+([1-9]+[0-9]*)')
    gse_list = ['GSE' + uid_regex.match(uid).group(1) for uid in result['IdList']]
    return(gse_list)
    