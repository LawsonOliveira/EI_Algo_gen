import mathutils
import math
import numpy as np



class RotTable:
    """Represents the rotation table"""

    # 3 first values: 3 angle values
    # 3 last values: SD values
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
    nucleotidlist=["AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"]
    def __init__(self):
        self.__Rot_Table = {}
        for dinucleotide in RotTable.__ORIGINAL_ROT_TABLE:
            #iniciate the values randomly
            self.__Rot_Table[dinucleotide] = [np.random.normal(RotTable.__ORIGINAL_ROT_TABLE[dinucleotide][0],RotTable.__ORIGINAL_ROT_TABLE[dinucleotide][3]),np.random.normal(RotTable.__ORIGINAL_ROT_TABLE[dinucleotide][1],RotTable.__ORIGINAL_ROT_TABLE[dinucleotide][4]),RotTable.__ORIGINAL_ROT_TABLE[dinucleotide][2]]

    def mut(self,k):
        #change the values asociated to the self.__Rot_Table[nucleotide] with nucleotide being the k_th nucleotide in the original table
        dinucleotide=RotTable.nucleotidlist[k]
        self.__Rot_Table[dinucleotide]=[np.random.normal(RotTable.__ORIGINAL_ROT_TABLE[dinucleotide][0],RotTable.__ORIGINAL_ROT_TABLE[dinucleotide][3]),np.random.normal(RotTable.__ORIGINAL_ROT_TABLE[dinucleotide][1],RotTable.__ORIGINAL_ROT_TABLE[dinucleotide][4]),RotTable.__ORIGINAL_ROT_TABLE[dinucleotide][2]]
    
    ###################
    # WRITING METHODS #
    ###################

    ###################
    # READING METHODS #
    ###################

    def getTwist(self, dinucleotide):
        return self.__Rot_Table[dinucleotide][0]

    def getWedge(self, dinucleotide):
        return self.__Rot_Table[dinucleotide][1]

    def getDirection(self, dinucleotide):
        return self.__Rot_Table[dinucleotide][2]

    ###################


print("on est passée ")



