'''
MODULE NAME:
    Dexs module

DESCRIPTION:
    Modulo que contiene la clase dexs

VERSION:
    1.0
    
AUTHORS
    Delgado Gutierrez Diana 
    Lopez Angeles Brenda E.
    Plascencia Perez Victor Ulises
    Ávila Silva Rogelio Lael
        
'''
import seaborn as sns
from pandas import *
import matplotlib.pyplot as plt

class dexs():
    def __init__(self, dexs_table):
        self.Ids = dexs_table.loc[:,['ID']]
        self.dexs_clustermap= sns.clustermap(dexs_table.loc[:,['control','treated']])
        plt.show()