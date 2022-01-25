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
    "TT": [35.62, 7.2,  154, ]
}  # Exemple de représentation
"""  Use a  tree to search in protein folding
We have 2 class: a tree and nodes
Rn: a node have a value, a successor function
Tree can chose his successor+ """
ÒÒÒÒÒÒÒÒ
# Some question: How to store the nodes values?

# Store data with a dictionnary of admissible interval

"""
class nodes(tree):
    # Set of data
    value = 0

    def init__(T):  # T is the ???  which represent all the possible interval
        nodes.rep = T  # all possible interval
        nodes.value = 0
        super().fit_node(self)

    def successor(self):
        # Renvoie la liste des successeurs
        res = []  # list of succesor
        for nuc in self.rep:

            for k in range(len(self.rep(nuc))):
                nextrep1 = getrep(self).copy()
                nextrep2 = getrep(self).copy()

                [a, b] = nextrep1[nuc][k]
                nextrep1[nuc][k] = [a, (a+b)/2]
                nextrep2[nuc][k] = [(a+b)/2, b]
        return res  # Non optimal: It creates way too much dictionnary

    def getrep(self):  # Renvoie le tableau de représentant associé
        return self.rep

"""


class node(object):  # Generic tree node
    "Generic tree node."

    def __init__(self):

        self.valeur = 0  # Valeur associé
        self.children = []  # Liste des children

    def __repr__(self):
        return self.name

    def add_child(self, node):
        assert isinstance(node, Tree)
        self.children.append(node)

    def getvalue(self):
        return self.valeur


t = node()
print(t.getvalue())


"""
   def fit_node(self):
        return
        # Fit a node

    def chosenode(self):
        pass  # Chose a node between all nodes possibles

    #  On a 4 étapes selection , expansion, simulaton, backpropagation


# Pour le tirage il faudrait prendre n valeurs au hasard  dans les intervalles et ensuite backpropager

# Pour """
