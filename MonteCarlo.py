from Traj3D import *
from random import *
import numpy as np
#uniform(-.5, .5)


class node:  # Generic tree node

    #"Generic tree node."
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
    nucleotidlist = ["AA", "AC", "AG", "AT", "CA", "CC", "CG",
                     "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]
    __ORIGINAL_INTERVALS = {}
    for key in __ORIGINAL_ROT_TABLE:
        __ORIGINAL_INTERVALS[key] = [
            [__ORIGINAL_ROT_TABLE[key][0]-np.sqrt(3)*__ORIGINAL_ROT_TABLE[key][3],
             __ORIGINAL_ROT_TABLE[key][0]+np.sqrt(3)*__ORIGINAL_ROT_TABLE[key][3]],
            [__ORIGINAL_ROT_TABLE[key][1]-np.sqrt(3)*__ORIGINAL_ROT_TABLE[key][4],
             __ORIGINAL_ROT_TABLE[key][1]+np.sqrt(3)*__ORIGINAL_ROT_TABLE[key][4]],
            [__ORIGINAL_ROT_TABLE[key][2]-np.sqrt(3)*__ORIGINAL_ROT_TABLE[key][5],
             __ORIGINAL_ROT_TABLE[key][2]+np.sqrt(3)*__ORIGINAL_ROT_TABLE[key][5]]
        ]
    # Exemple de représentation
    # """  Use a  tree to search in protein folding
    # We have 2 class: a tree and nodes
    # Rn: a node have a value, a successor function
    # Tree can chose his successor+ """
    # Some question: How to store the nodes values?

    # Store data with a dictionnary of admissible interval

    # """
    # class nodes(tree):
    #     # Set of data
    #     value = 0

    #     def init__(T):  # T is the ???  which represent all the possible interval
    #         nodes.rep = T  # all possible interval
    #         nodes.value = 0
    #         super().fit_node(self)

    #     def successor(self):
    #         # Renvoie la liste des successeurs
    #         res = []  # list of succesor
    #         for nuc in self.rep:

    #             for k in range(len(self.rep(nuc))):
    #                 nextrep1 = getrep(self).copy()
    #                 nextrep2 = getrep(self).copy()

    #                 [a, b] = nextrep1[nuc][k]
    #                 nextrep1[nuc][k] = [a, (a+b)/2]
    #                 nextrep2[nuc][k] = [(a+b)/2, b]
    #         return res  # Non optimal: It creates way too much dictionnary

    #     def getrep(self):  # Renvoie le tableau de représentant associé
    #         return self.rep

    # """

    def __init__(self, table=__ORIGINAL_ROT_TABLE, interval=__ORIGINAL_INTERVALS):

        self.__Rot_Table = table  # to complete

        self.__valeur = 0  # score of the function
        self.__n = 0          # number of time we chose this node
        self.__h = 0          # Height of the tree
        self.__Childs = []
        self.__intervals = node.__ORIGINAL_INTERVALS

    def add_child(self, node):
        self.__Childs.append(node)

    def writeRot_Table(self, dict):
        # dict is the dictionary that will inserted in the Rot_Table of self
        # has the shape {"AA": [154, 7, -154], "GC":[40, 5,  180, ]}
        for key in dict:
            assert key in self.__Rot_Table
            self.__Rot_Table[key] = dict[key]

    def calculateD(self, seq):  # Calculate Traj 3D
        traj = Traj3D()
        traj.compute(seq, self.__Rot_Table)
        xyz = np.array(traj.getTraj())
        x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
        return np.sqrt(x[-1]**2 + y[-1]**2 + z[-1]**2)

    def getvalue(self):  # return the score of the function
        return self.__valeur

    def writeValeur(self, idk):  # don't know what to put here
        self.__valeur = idk

    def getTable(self):
        return self.__Rot_Table

    def getn(self):
        return self.__n

    def actualizen(self, k):
        self.__n = k

    def geth(self):
        return self.__h

    def actualizeh(self, m):
        self.__h = m

    def getinterval(self):
        return self.__intervals

    def getintervalspec(self, nuc, angle):
        return self.__intervals[nuc][angle]

    def actualiseinterval(self, nuc, angle, value):
        self.__intervals[nuc][angle] = value

    def getoriginalrotable(self):
        return self.__ORIGINAL_ROT_TABLE

    def getchild(self):
        return self.__Childs


# """
#    def fit_node(self):
#         return
#         # Fit a node
#     def chosenode(self):
#         pass  # Chose a node between all nodes possibles
#     #  On a 4 étapes selection , expansion, simulaton, backpropagation
# # Pour le tirage il faudrait prendre n valeurs au hasard  dans les intervalles et ensuite backpropager
# # Pour """
