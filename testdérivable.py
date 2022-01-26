import argparse
from RotTable import *
from Traj3D import *
# La geule de la matrice est pas clair donc inexploitable en fait....

from math import sqrt

# hyperparametre
pasderive = 10**-3


def calculderiveepartielle(f, L, pasderive, x, epsilon=1):
    # X est l'indice de la  variable , le pas est trivialement  le pas et f la fonction à n variable
    """Renvoie la valeur de la derivée partielle en x  """
    # Ici, pour diminuer la complexité, on calcule la derive partielle in place, pour cela on ajoute directement à la liste L
    # epsilon est le caractère positif ou négatif de la dérivée partielle
    a1 = f(L)
    for element in L:
        for k in range(3):

            L.__Rot_Table[k] = L.__Rot_Table[k]+epsilon*pasderive
            aapres = f(L)  # c'est le f(x+h)
            L.__Rot_Table[k] = L.__Rot_Table[k]+epsilon*pasderive

    return (aapres-a1)/pasderive


def norme(L):
    """ Renvoie la norme euclidienne des n uplets (a,b,c,d): a²+b²+c²+d² ...
    """
    resultat = 0
    for k in L:
        resultat += k**2
    return resultat


def normematrix(L):
    resultat = 0
    for l in L:
        for k in l:
            resultat += k**2
    return resultat


def gradientpos(f, L, pasderive, epsilon=1):
    """On renvoie le gradient des dérivées partielles dans l'ordre usuelle"""
    # L désigne la liste qui contient les n variables de f, """
    resultat = []  # c'est le gradient
    for indicevariable in range(len(L)):
        resultat.append(calculderiveepartielle(f, L, pasderive,
                        indicevariable, epsilon))   # on ajouter df/dx
    return resultat


def gradienttest(f, L, pasderive, epsilon=1):
    resultat = []

    for element in L._RotTable__Rot_Table:
        for k in range(3):
            a1 = f(L)

            L._RotTable__Rot_Table[element][k] = L._RotTable__Rot_Table[element][k] + \
                epsilon*pasderive
            aapres = f(L)  # c'est le f(x+h)
            L._RotTable__Rot_Table[element][k] = L._RotTable__Rot_Table[element][k] - \
                epsilon*pasderive
            res = []
            for k in range(len(aapres)):
                res.append((aapres[k]-a1[k]/pasderive))
            resultat.append(res)

    return resultat


# ici on vérifie si le gradient est positif
parser = argparse.ArgumentParser()
parser.add_argument("--filename", help="input filename of DNA sequence")
parser.parse_args()
args = parser.parse_args()


def main():

    a = RotTable()
    traj = Traj3D()

    def ftraj(rot_table):
        traj.compute(
            "AAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCCAGTAAACGAAAAAACCGCCTGGGGAGGCGGTTTAGTCGAAGGTTAAGTCAG", rot_table)
        return traj.getmatrix()

    f = ftraj

    # print(f(a))
    n = 16  # Nombre de element possible
    intadm = {}  # Ensemble des intervalles admissibles
    """  print(a._RotTable__Rot_Table)
    for element in a._RotTable__Rot_Table:
        intadm[element] = {}
        for k in range(3):

            (ak, bk) = (intadm[k]+intadm[k+3], intadm[k]-intadm[k+3])
            intadm[element][k] = (ak, bk)

    # On a ainsi une liste des interval"""
    pasderive = 10**-3
    sample = 50
    test = []
    epsilon = 10**-2
    marqueur = True
    for k in range(sample):
        a = RotTable()

        aplus = gradienttest(f, a, pasderive, 1)
        amoins = gradienttest(f, a, pasderive, -1)

        averif = []
        for k in range(len(aplus)):
            ligne = []
            plusa = np.array(aplus[k])
            amoinsa = np.array(amoins[k])
            for i in range(len(plusa[k])):

                ligne.append(plusa[0][i]-amoinsa[0][i])

            averif.append(ligne)
        """print("averif", averif)
        print(averif, "de norme:", normematrix(averif))
        print(" et LA DERIVE PARTIELLE EST BIEN","""
        marqueur = marqueur and normematrix(averif) < epsilon
        if k//100 == 0:
            print(k)
            print(marqueur)

    print(marqueur, " ET LE RESULTAT EST ")

    inttest = []


if __name__ == "__main__":
    marqueur = True
    main()
