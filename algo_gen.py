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

    ###################


class population():

    def __init__(self, selection, mutation, fit):
        self.pop = []
        self.selection = selection
        self.mutation = mutation
        self.fit = fit


def initialisation_population():  # Initialize the population ,
    pass


critere = True  # a modifier


def Algogenetique():
    # ici
    pop = initialisation_population()

    while critere:  # Critere has to be defined
        va = pop.evaluation()

        pop1 = pop.selection(pop, va)
        pop1 = pop.mutation()
        pop = pop1

    return pop
