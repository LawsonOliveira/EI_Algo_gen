# CARE WE ARE MAYBE NOT WORKING WITH ROT TABLE
import numpy as np
from math import sqrt
from RotTable import RotTable
from numpy import angle
# from EI_Algo_gen.RotTable import RotTable
from MonteCarlo import *
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
# Hyperparameter

seq = "AAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCCAGTAAACGAAAAAACCGCCTGGGGAGGCGGTTTAGTCGAAGGTTAAGTCAGTTGGGGACTGCTTAACCGGGTAACTGGCTTGGTGGAGCACAGATACCAAATACTGTCCTTCTAGTGTAGCCGCAGTTAGGCCACCACTTCAAGAACTCTTAATATCTCAATCCACCTTGTCCAGTTACCAGTGGCTGCTGCCAGTGGCGCTTTGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAGCCGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTATGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAACGGCGCAGCCTTTTCCTGGTTCTCGTTTTTTGCTCACATGTTTCTTTTGGCGTTATCCCCTGATTCTGTGGATAACCGCATCTCCGCTTTTGAGTGAGCAGACACCGCTCGCCGCAGCCGAACGACCGAGTGTAGCGAGTCAGTGAGCGAGGAAGCGGAAGAGCGCCGGAACGTGCATTTTCTCCTTACGCATCTGTGCGGCATTTCACATCGGACATGGTGCGCTTTCCATACAATTCGTACTGATGCCGCATAGTTAAGCCAGTATACACTCCGCTATCGCTACGTGACTGGTTCAGGGCTTCGCCCCGAAACCCCCTGACGCGCCCTGAGGGGCTTGTCTGCTCCCGGCATCCGCTCACAGACAAGCTGTTACCGTCTCCGGGAGCTGTATGTGTCAGAGGTTTTCACCGTCATCCCCGAAGCGTGCGA"
nbsamplesmallh = 100  # Parameter for sample for the first depth
nbsamplehugeh = 50  # first parameter for the sample for huge depth
barrier = 16  # How much we have huge or bad sample
# Three Main part of Monte Carlo : selection, expansion, backpropagation


def fit(dict, seq=seq):
    a = RotTable()
    a.newTable(dict)

    traj = Traj3D()
    traj.compute(seq, a)
    xyz = np.array(traj.getTraj())
    x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    # BEWARE THE MINUS WE  HAVE TO MAXIMIZE THE FUNCTION
    return -np.sqrt(x[-1]**2 + y[-1]**2 + z[-1]**2)


def selection(node, K=1):
    childlist = node.getChilds()  # list of children
    Nchild = len(childlist)  # The len of the child to have
    N = node.getn()  # THe value of the father

    def valuebandit(child):  # Return the Bandit value
        score = child.getvalue()  # Value do to the fit function
        nbseen = child.getn()+1  # Nb of this time the son has been visited

        # ATTENTION ON RECUPERE UN SCORE MAXIMAL ICI ET PAS UN SCORE MINIMAL, LA FORMULE DAN NOTRE CAS EST FAUSSE
        return score/nbseen + K*np.sqrt(3/2*np.log(N)/nbseen)
        # Constnate K est louche
    # return the best in childlist who has the best value over bandit
    childchoisi = max(childlist, key=valuebandit)
    # Don't forget to add the time we add the chose of child
    childchoisi.writen(childchoisi.getn()+1)
    return childchoisi


def sampleininterval(sample, interval):  # Check if a sample is in the interval or not
    marq = True
    for nuc in sample:
        for value in range(3):
            if interval[nuc][value][0] > sample[nuc][value] or interval[nuc][value][1] < sample[nuc][value]:
                return False

    return True


# m mean how we gonna split the interval in many parts, we create here m Child
def createchild(node1, m=10):
    #  On créer m enfants
    nbnucleotide = 16  # nb of nucleotide
    nbangle = 3
    h = node1.geth() + 1  # The depth of the child
    # Here two things we want:
    hnew = h % (nbnucleotide*nbangle)

    nuc = node.nucleotidlist[hnew//3]  # The nucleotide we gonna have new

    # Now we have to create the childe
    # __intervals[key]

    node1.actualizen(node1.getn()+m)
    anglestudied = hnew % 3
    b = node1.getinterval()[nuc][anglestudied][1]
    a = node1.getinterval()[nuc][anglestudied][0]
    marqueur = True  # Marqueur checks  if sample is in the interval on not
    for i in range(m):
        n_nodes = node()
        n_nodes.actualizeh(h)
        # we copy the dictionnary we are looking at , # The best would be to have a list, and we do it directy on the dictionnary ...
        n_nodes.actualiseinterval(
            nuc, anglestudied, [a + (b-a)*i/m, a + (b-a)*(i+1)/m])

        if h < barrier:
            evaluate(n_nodes, nbsamplesmallh)
        if h > barrier:  # We modify, cuz the first one are pretty important

            evaluate(n_nodes, nbsamplehugeh)
        # To prevent, if the sample we got before are as good as before

        rottable = node1.getTable()

        interval = n_nodes.getinterval()
        if marqueur and sampleininterval(rottable, interval):
            marqueur = False  # We don't need to search any further
            value = node1.getvalue()
            if value > n_nodes.getvalue():
                n_nodes.writeValeur(value)
                n_nodes.writeRot_Table(rottable)
        # End of prevent,
        node1.add_child(n_nodes)  # We add finally the children


def evaluate(node, nbsample=30):  # evaluate a node, by taking a lot of children

    # Nb of sample we get
    min = (-1)
    Rot_table = {}
    h = node.geth()  # getH

    interval = node.getinterval()  # we get back the interval
    for nuc in interval:

        Rot_table[nuc] = {}

    for k in range(nbsample):  # Question is how much sample we want to study

        samplestudied = {}

        for nuc in interval:  # Pour chaque nucléotide
            samplestudied[nuc] = []

            for anglestudied in interval[nuc]:  # Pour chaque angle
                #print("anglestudied", anglestudied, "interval", nuc)
                lbound = anglestudied[0]
                hbound = anglestudied[1]
                assert lbound <= hbound, ' intervalle pas bien '
                ak = np.random.uniform(lbound, hbound)

                samplestudied[nuc].append(ak)

            # May be unusufull

            # Pour se conformer à la fit fonction
            samplestudied[nuc] = samplestudied[nuc] + \
                Rot_Table[nuc][4:]

        # Evaluation of the sample, SEE GENETIC ALGORITHM
        #print("samplestudied", samplestudied)
        value = fit(samplestudied, seq)

        if (min == -1) or value > min:
            node.writeValeur(value)  # We write the new value
            min = value
            # We modify the Rot_table right now
            node.writeRot_Table(samplestudied)


# NON FAIT CREATE CHILD, EVALUATE
def expansion(node):
    # here node doesn't have child
    assert node.getchild() == [], 'No a feather'
    # print(node)
    # Pour l'expansion , on ouvre k enfants , on les évalues, on récupère le maximum et le renvoie pour mettre à jour les enfants
    N = 10  # Number of child we create
    valevaluate = []  # list of score of each son
    h = node.geth()

    # We create a child , IT DEPENDS WHAT WE DO HERE
    createchild(node, N)  # Create child creates N child to the "node" node

    childlist = node.getchild()
    for child in childlist:
        ak = child.getvalue()
        valevaluate.append(ak)
    # Don't forget to add the child

    m = max(valevaluate)
    u = node.getvalue()
    if u == 0 or m > u:

        node.writeValeur(m)   # We modify the value of the node
    return max(u, m)  # We return the value of the node


def backpropagation(node):  # Algorithme  finale de backpropagation
    visiter = node
    securityrope = [visiter]  # To hinder recursivity

    childvisiter = visiter.getchild()
    while childvisiter != []:
        # It has children

        nvavister = selection(visiter)  # We chose a new children

        visiter = nvavister  # Then the new node we are looking is visiter
        childvisiter = nvavister.getchild()
        securityrope.append(visiter)
    term = securityrope[-1]  # The last one is a feather
    value = expansion(term)  # EXPANSIONNNN
    # Time to BACKPEDAL
    for node in securityrope:
        nodeval = node.getvalue()  # we take the value of the node we are looking at

        if nodeval > value:
            node.__value = value
    return  # We return nothing, we don't care about it


def compute(root, nbit, critere=10**-3):
    # Root: C.I
    """NBIT: NOMBRE ITERATION max 
    CRITERE POUR SAVOIR SI ON EST PROCHE DE LA DISTANCE OU NON """

    value = critere+1  # Value of the first root, MUST NOT BE 0 WATCH OUT
    nbiteration = 0
    print("initialisation", value)
    print()
    while nbiteration < nbit and abs(value) > critere:
        backpropagation(root)  # We backpropagation like always
        value = root.getvalue()
        print(" nbiteration", nbiteration)
        print("value", value)
        nbiteration += 1

    # on a donc le graphe final

    #  On récupère la solution
    def recupsol(node):
        a = node
        while a.getchild() != []:
            childliste = a.getchild()
            nodee = a
            max = -1
            for child in childliste:
                if max == -1 or child.getvalue() > max:
                    max = child.getvalue()
                    nodee = child

            # Here we have the child with the best value
            a = nodee
        return a.getTable()

    dive = root
    print(" le plus dur est passée, on récupère le mnimum ")
    return recupsol(dive)  # Return the best sample we ever had


# APPEND  copie superficielle


def main():
    noeud = node()
    nbit = 10  # HYPERPARAMETRE
    best = compute(noeud, nbit)
    #print(noeud.getvalue(), "his value")
    print(noeud.getvalue(), "la value")
    print(noeud.getchild(), "les momes")

    traj = Traj3D()
    print(best)  # BEST IS AN ARRAY
    a = Rot_Table()
    a.newTable(best)
    traj.compute(seq, a)

    # traj.compute(traj,seq,best)
    traj.draw("MONRESULTAT.png")
    print("pitié ça marche")


if __name__ == "__main__":

    a = node()
    Rot_Table = a.getoriginalrotable()  # TAKE THIS AS GLOBAL VARIABLE
    main()
