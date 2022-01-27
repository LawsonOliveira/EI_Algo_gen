import numpy as np
import random

class Chromosome:
    __ORIGINAL_ROT_TABLE = {
        "AA": [35.62, 7.2, -154,      0.06,  0.6, 0],
        "AC": [34.4, 1.1,  143,      1.3,  5, 0],
        "AG": [27.7, 8.4,    2,      1.5,  3, 0],
        "AT": [31.5, 2.6,    0,      1.1,  2, 0],
        "CA": [34.5, 3.5,  -64,      0.9, 34, 0],
        "CC": [33.67, 2.1,  -57,      0.07,  2.1, 0],
        "CG": [29.8, 6.7,    0,      1.1,  1.5, 0],
        "CT": [27.7, 8.4,   -2,      1.5,  3, 0],
        "GA": [36.9, 5.3,  120,      0.9,  6, 0],
        "GC": [40, 5,  180,      1.2,  1.275, 0],
        "GG": [33.67, 2.1,   57,      0.07,  2.1, 0],
        "GT": [34.4, 1.1, -143,      1.3,  5, 0],
        "TA": [36, 0.9,    0,      1.1,  2, 0],
        "TC": [36.9, 5.3, -120,      0.9,  6, 0],
        "TG": [34.5, 3.5,   64,      0.9, 34, 0],
        "TT": [35.62, 7.2,  154,      0.06,  0.6, 0]
    }

    def __init__(self,genes=None):
        if genes==None:
            self.__genes = {}
            for dinucleotide in self.__ORIGINAL_ROT_TABLE:
                #iniciate the values randomly
                self.__genes[dinucleotide] = [np.random.normal(self.__ORIGINAL_ROT_TABLE[dinucleotide][0],self.__ORIGINAL_ROT_TABLE[dinucleotide][3]),np.random.normal(self.__ORIGINAL_ROT_TABLE[dinucleotide][1],self.__ORIGINAL_ROT_TABLE[dinucleotide][4]),self.__ORIGINAL_ROT_TABLE[dinucleotide][2]]
                #self.__genes[dinucleotide] = [np.random.uniform(self.__ORIGINAL_ROT_TABLE[dinucleotide][0]-self.__ORIGINAL_ROT_TABLE[dinucleotide][3],self.__ORIGINAL_ROT_TABLE[dinucleotide][0]+self.__ORIGINAL_ROT_TABLE[dinucleotide][3]),np.random.uniform(self.__ORIGINAL_ROT_TABLE[dinucleotide][1]-self.__ORIGINAL_ROT_TABLE[dinucleotide][4],self.__ORIGINAL_ROT_TABLE[dinucleotide][1]+self.__ORIGINAL_ROT_TABLE[dinucleotide][4]),np.random.uniform(self.__ORIGINAL_ROT_TABLE[dinucleotide][2]-self.__ORIGINAL_ROT_TABLE[dinucleotide][5],self.__ORIGINAL_ROT_TABLE[dinucleotide][2]+self.__ORIGINAL_ROT_TABLE[dinucleotide][5])]
        else:
            self.__genes=genes
        self.__nucleotidlist=["AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"]
        self.__nb_genes = len(self.__ORIGINAL_ROT_TABLE)
    
    def __str__(self):
        return str(self.__genes)

    def crossover(self,chromosome1, chromosome2):
        #we initialize chr_os1 and chr_os2: chromosome offsprings
        enfant1={}
        enfant2={}
        for dinucleotide in self.__nucleotidlist:
            enfant1[dinucleotide]=random.choice([chromosome1.get_gene(dinucleotide),chromosome2.get_gene(dinucleotide)])
            enfant2[dinucleotide]=random.choice([chromosome1.get_gene(dinucleotide),chromosome2.get_gene(dinucleotide)])
        return enfant1,enfant2

    def apply_mutation(self):
        dinucleotide=random.randint(0,len(self.__nucleotidlist)-1)
        dinucleotide=self.__nucleotidlist[dinucleotide]
        self.__genes[dinucleotide] = [np.random.uniform(self.__ORIGINAL_ROT_TABLE[dinucleotide][0]-self.__ORIGINAL_ROT_TABLE[dinucleotide][3],self.__ORIGINAL_ROT_TABLE[dinucleotide][0]+self.__ORIGINAL_ROT_TABLE[dinucleotide][3]),np.random.uniform(self.__ORIGINAL_ROT_TABLE[dinucleotide][1]-self.__ORIGINAL_ROT_TABLE[dinucleotide][4],self.__ORIGINAL_ROT_TABLE[dinucleotide][1]+self.__ORIGINAL_ROT_TABLE[dinucleotide][4]),np.random.uniform(self.__ORIGINAL_ROT_TABLE[dinucleotide][2]-self.__ORIGINAL_ROT_TABLE[dinucleotide][5],self.__ORIGINAL_ROT_TABLE[dinucleotide][2]+self.__ORIGINAL_ROT_TABLE[dinucleotide][5])]

    def upd_chr(self,genes):
        self.__genes = genes    
    
    def get_size(self):
        return len(self.__genes)
    
    def get_chr(self):
        return self.__genes
    
    def get_gene(self, dinucleotide):
        return self.__genes[dinucleotide]

    def getTwist(self, dinucleotide):
        return self.__genes[dinucleotide][0]

    def getWedge(self, dinucleotide):
        return self.__genes[dinucleotide][1]

    def getDirection(self, dinucleotide):
        return self.__genes[dinucleotide][2]

