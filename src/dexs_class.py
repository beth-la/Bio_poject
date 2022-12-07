import seaborn as sns
from pandas import *
import matplotlib.pyplot as plt

class dexs():
    def __init__(self, dexs_table):
        self.dexs_Ids = dexs_table.loc[:,['ID']]
        self.dexs_cluster_map= sns.clustermap(dexs_table.loc[:,['control','treated']])
        plt.show()
