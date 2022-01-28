import sys  
from pathlib import Path  
file = Path(__file__).resolve()  
package_root_directory = file.parents[1]  
sys.path.append(str(package_root_directory))

import numpy as np
import random

# This class represents each chromosome
class Chromosome:
    # Original rotation table
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

    def __init__(self,chr=None):
        # Chr is a collection of gene of a chromosome
        if chr==None:
            self.__chr = {}
            for dinucleotide in self.__ORIGINAL_ROT_TABLE:
                # Iniciate the values randomly with a uniform distribution 
                self.__chr[dinucleotide] = [np.random.uniform(self.__ORIGINAL_ROT_TABLE[dinucleotide][0]-self.__ORIGINAL_ROT_TABLE[dinucleotide][3],self.__ORIGINAL_ROT_TABLE[dinucleotide][0]+self.__ORIGINAL_ROT_TABLE[dinucleotide][3]),np.random.uniform(self.__ORIGINAL_ROT_TABLE[dinucleotide][1]-self.__ORIGINAL_ROT_TABLE[dinucleotide][4],self.__ORIGINAL_ROT_TABLE[dinucleotide][1]+self.__ORIGINAL_ROT_TABLE[dinucleotide][4]),np.random.uniform(self.__ORIGINAL_ROT_TABLE[dinucleotide][2]-self.__ORIGINAL_ROT_TABLE[dinucleotide][5],self.__ORIGINAL_ROT_TABLE[dinucleotide][2]+self.__ORIGINAL_ROT_TABLE[dinucleotide][5])]
        else:
            # Creates a new chromosome with the given chr
            self.__chr=chr
        self.__nucleotidlist=["AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"]
    
    def __str__(self): 
        # Sends a chromosome string
        return str(self.__chr)

    def crossover(self,chromosome1, chromosome2):
        # We use couple chromosome1 and chromosome2 to create two children (news chromosomes)
        enfant1={}
        enfant2={}
        for dinucleotide in self.__nucleotidlist:
            # Each gene of each child is chosen randomly between the gene of the chromosome1 and the gene of the chromosome2
            shuffle=[chromosome1.get_gene(dinucleotide),chromosome2.get_gene(dinucleotide)]
            random.shuffle(shuffle)
            enfant1[dinucleotide]=shuffle[0]
            enfant2[dinucleotide]=shuffle[1]
        return enfant1,enfant2

    def apply_mutation(self):
        # Mutates the value of self.__chr with a uniform distribution 
        dinucleotide=random.randint(0,len(self.__nucleotidlist)-1)
        dinucleotide=self.__nucleotidlist[dinucleotide]
        self.__chr[dinucleotide] = [np.random.uniform(self.__ORIGINAL_ROT_TABLE[dinucleotide][0]-self.__ORIGINAL_ROT_TABLE[dinucleotide][3],self.__ORIGINAL_ROT_TABLE[dinucleotide][0]+self.__ORIGINAL_ROT_TABLE[dinucleotide][3]),np.random.uniform(self.__ORIGINAL_ROT_TABLE[dinucleotide][1]-self.__ORIGINAL_ROT_TABLE[dinucleotide][4],self.__ORIGINAL_ROT_TABLE[dinucleotide][1]+self.__ORIGINAL_ROT_TABLE[dinucleotide][4]),np.random.uniform(self.__ORIGINAL_ROT_TABLE[dinucleotide][2]-self.__ORIGINAL_ROT_TABLE[dinucleotide][5],self.__ORIGINAL_ROT_TABLE[dinucleotide][2]+self.__ORIGINAL_ROT_TABLE[dinucleotide][5])]

    def upd_chr(self,chr):
        # Updates the chromosome
        self.__chr = chr    
    
    def get_size(self):
        # Return the size of the chromosome
        return len(self.__chr)
    
    def get_chr(self):
        # Return the chromosome
        return self.__chr
    
    def get_gene(self, dinucleotide):
        # Return the gene in the dinucleotide position of the chromosome
        return self.__chr[dinucleotide]

    def getTwist(self, dinucleotide):
        # Return the twist of the dinucleotide in the chromosome
        return self.__chr[dinucleotide][0]

    def getWedge(self, dinucleotide):
        # Return the wedge of the dinucleotide in the chromosome
        return self.__chr[dinucleotide][1]

    def getDirection(self, dinucleotide):
        # Return the direction of the dinucleotide in the chromosome
        return self.__chr[dinucleotide][2]

