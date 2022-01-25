from math import sqrt
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


# Three Main part of Monte Carlo : selection, expansion, backpropagation


def selection(node, K=1):
    childlist = node.getchildren()  # list of children
    Nchild = len(childlist)  # The len of the child to have
    N = node.score()  # THe value of the father

    def valuebandit(child):  # Return the Bandit value
        score = child.score()  # Value do to the fit function
        nbseen = child.getn()  # Nb of this time the son has been visited

        # ATTENTION ON RECUPERE UN SCORE MAXIMAL ICI ET PAS UN SCORE MINIMAL, LA FORMULE DAN NOTRE CAS EST FAUSSE
        return score/nbseen + K*sqrt(3/2*ln(N)/nbseen)
        # Constnate K est louche
    # return the best in childlist who has the best value over bandit
    childchoisi = max(childlist, key=valuebandit)
    childchoisi.modn = nbseen+1  # Don't forget to add the time we add the chose of child
    return childchoisi


# NON FAIT CREATE CHILD, EVALUATE
def expansion(node):
    # here node doesn't have child
    assert node.child == []
    # Pour l'expansion , on ouvre k enfants , on les évalues, on récupère le maximum et le renvoie pour mettre à jour les enfants
    N = 100  # Number of child we create
    valevaluate = []  # list of score of each son
    for k in range(N):
        # We create a child , IT DEPENDS WHAT WE DO HERE
        child = createchild(node)

        ak = evaluate(child)

        child.score = ak  # The value we got
        child.depth = node.depth+1  # Add the depth

        valevaluate.append(ak)

        # Don't forget to add the child
        node.ajouterchild(child)
    m = max(valevaluate)
    node.modscore = m  # We modify the value of the node
    return m  # We return the value of the node
    # Here, the question is how we expand ?


def backpropagation(node):  # Algorithme  finale de backpropagation
    securityrope = []  # To hinder recursivity
    visiter = node

    while visiter.child != []:
        securityrope.append(visiter)  # It has children

        nvavister = selection(visiter)  # We chose a new children
        visiter = nvisiter  # Then the new node we are looking is visiter

    term = securityrope[-1]  # The last one is a feather
    value = expansion(node)  # EXPANSIONNNN
    # Time to BACKPEDAL
    for node in securityrope:
        nodeval = node.getvalue()  # we take the value of the node we are looking at

        if nodeval > value:
            node.modvalue = value
    return  # We return nothing, we don't care about it


def compute(root, nbit, critere=10**-3):
    # Root: C.I
    """NBIT: NOMBRE ITERATION max 
    CRITERE POUR SAVOIR SI ON EST PROCHE DE LA DISTANCE OU NON """

    value = root.getvalue()  # Value of the first root, MUST NOT BE 0 WATCH OUT
    nbiteration = 0

    while nbiteration < nbit and value > critere:
        backpropagation(root)  # We backpropagation like always

    # on a donc le graphe final

    #  On récupère la solution

    dive = root

    while dive.child[] != []:

        # CARE MIN OR MAX , evaluer fonction non fait
        dive = max(dive.child, key=evaluerfonction)

    return dive.sampleevalue()  # Return the best sample we ever had


# APPEND  copie superficielle
