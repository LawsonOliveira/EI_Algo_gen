from os import popen
from re import A


nbseqnucleotique = 16  # Nombre de séquence possible
nbangle = 3  # Nombre d'angle


class chromosome:
    pop = 0
    # Valeur  absolue des angles admissibles
    maxanglevalue = {
        "AA": [0.06,  0.6, 0],
        "AC": [1.3,  5, 0],
        "AG": [1.5,  3, 0],
        "AT": [1.1,  2, 0],
        "CA": [0.9, 34, 0],
        "CC": [0.07,  2.1, 0],
        "CG": [1.1,  1.5, 0],
        "CT": [1.5,  3, 0],
        "GA": [0.9,  6, 0],
        "GC": [1.2,  1.275, 0],
        "GG": [0.07,  2.1, 0],
        "GT": [1.3,  5, 0],
        "TA": [1.1,  2, 0],
        "TC": [0.9,  6, 0],
        "TG": [0.9, 34, 0],
        "TT": [0.06,  0.6, 0]
    }

    def __init__(self):
        #  Initialize normal value
        self.__Rot_Table = {
            "AA": [35.62, 7.2, -154],
            "AC": [34.4, 1.1,  143],
            "AG": [27.7, 8.4,    2],
            "AT": [31.5, 2.6,    0],
            "CA": [34.5, 3.5,  -64],
            "CC": [33.67, 2.1,  -57],
            "CG": [29.8, 6.7,    0],
            "CT": [27.7, 8.4,   -2],
            "GA": [36.9, 5.3,  120],
            "GC": [40, 5,  180],
            "GG": [33.67, 2.1,   57],
            "GT": [34.4, 1.1, -143],
            "TA": [36, 0.9,    0],
            "TC": [36.9, 5.3, -120],
            "TG": [34.5, 3.5,   64],
            "TT": [35.62, 7.2,  154]
        }

    # Fonction à avoir : récupérer
    def getTwist(self, dinucleotide):
        return self.__Rot_Table[dinucleotide][0]

    def getWedge(self, dinucleotide):
        return self.__Rot_Table[dinucleotide][1]

    def getDirection(self, dinucleotide):
        return self.__Rot_Table[dinucleotide][2]

    def traj3D.compute(self, dna_seq, self.rot_table)
    ###################


class population():

    def __init__(self, selection, mutation, fit):
        self.pop = []
        self.selection = selection
        self.mutation = mutation
        self.fit = fit


def initialisation_population():  # Initialize the population ,
    pass


class GA:
    """
    Classe permettant d'éxécuter un algorithme génétique sur un problème d'optimisation.
    """

    def __init__(self, pop_cap, nb_var, fitness_fnct):
        """
        Sauvegarde les hyper-paramètres de l'algorithme, et crée la population initiale.

        Paramètres :
        - pop_cap (int) : Nombre d'individus de la population
        - nb_var (int) : Taille du génome des individus
        - fitness_fnct (function) : Fonction objective à minimiser
        """
        # Sauvegarde des hyper-paramètres
        self.pop_cap = pop_cap
        self.nb_var = nb_var
        self.fitness_fnct = fitness_fnct

    def ini_pop(self):
        pass

    def selection(self, ind1, ind2):
        pass

    def mutation(self, ind, mut_rate,
                 bound_inf=None, bound_sup=None):

        return ind

    def evaluation(self):
        """ évalue le génome actuelle"""

    def do_gen(self):
        """
        Calcule la génération N+1 à partir de la génération N.
        """
        pass

    def compute():
        pop = initialisation_population()
        while critere:  # Critere has to be defined
            va = pop.evaluation()

            pop1 = pop.selection(pop, va)
            pop1 = pop.mutation()
            pop = pop1

        return pop

# Liste à trois dimensions


def fit1(trajchromosome1, trajchromosome2):
    def distance(point1, point2):
        return (point1[0]-point2[0])**2 + (point1[1]-point2[1])**2 + (point1[2]-point2[2])**2
    # Distance associé aux normes euclidiennes
    return distance(trajchromosome1, trajchromosome2)


# Méta donnée Pc
def selection(pop):
    # Liste de population
    pop = sorted(words, key=fit)  # On ordonne selon la fonction len
    evaluerpop(pop, fit)
