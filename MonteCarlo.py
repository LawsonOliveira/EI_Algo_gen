from Traj3D import *
from random import *

class node:  # Generic tree node

    #"Generic tree node."
    __ORIGINAL_ROT_TABLE = {
    "AA": [35.62, 7.2, -154, ],
    "AC": [34.4, 1.1,  143, ],
    "AG": [27.7, 8.4,    2, ],
    "AT": [31.5, 2.6,    0, ],
    "CA": [34.5, 3.5,  -64, ],
    "CC": [33.67, 2.1,  -57, ],
    "CG": [29.8, 6.7,    0, ],
    "CT": [27.7, 8.4,   -2, ],
    "GA": [36.9, 5.3,  120, ],
    "GC": [40, 5,  180, ],
    "GG": [33.67, 2.1,   57, ],
    "GT": [34.4, 1.1, -143, ],
    "TA": [36, 0.9,    0, ],
    "TC": [36.9, 5.3, -120, ],
    "TG": [34.5, 3.5,   64, ],
    "TT": [35.62, 7.2,  154, ]}  
    
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

    def __init__(self):

        self.__Rot_Table={} #to complete

        
        self.__valeur=0     #to complete
        self.__n=0          
        self.__N=0
        self.__Childs=[]

    def add_child(self, node):
        self.__Childs.append(node)

    def writeRot_Table(self,dict):
        #dict is the dictionary that will inserted in the Rot_Table of self
        #has the shape {"AA": [154, 7, -154], "GC":[40, 5,  180, ]}
        for key in dict:
            self.__Rot_Table[key]=dict[key]

    def writeValeur(self,idk): #don't know what to put here
        pass
    
    def calculateD(self,seq):
        traj=Traj3D()
        traj.compute(seq,self.__Rot_Table)        
        xyz = np.array(traj.getTraj())  
        x, y, z = xyz[:,0], xyz[:,1], xyz[:,2]
        return np.sqrt(x[-1]**2 + y[-1]**2 + z[-1]**2)

    def getvalue(self):
        return self.__valeur
    
    def actualizeN(self):
        pass


test= node()
print(test.getvalue())

# """
#    def fit_node(self):
#         return
#         # Fit a node
#     def chosenode(self):
#         pass  # Chose a node between all nodes possibles
#     #  On a 4 étapes selection , expansion, simulaton, backpropagation
# # Pour le tirage il faudrait prendre n valeurs au hasard  dans les intervalles et ensuite backpropager
# # Pour """