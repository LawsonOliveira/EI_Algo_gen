
def crossover2(bests):
    index={i:i for i in range(len(bests))}
    fil={}
    fils=[]
    taille=10*len(bests)
    while len(index)>1:
        for n in range(taille):
            i=random.choice(list(index.values()))
            index.pop(i)
            j=random.choice(list(index.values()))
            coupe=random.randint(0,len(bests)-1)
            counter=0
            for dinucleotide in bests[i].getTable():
                if counter<=coupe:
                    fil[dinucleotide]=bests[i].getTable()[dinucleotide]
                if counter>coupe:
                    fil[dinucleotide]=bests[j].getTable()[dinucleotide]
                counter+=1
            a=RotTable()
            a.newTable(fil)
            fils.append(a)
            fil={}
            index[i]=i
        index.pop(i)
        index.pop(j)
    return fils        


def crossover3(bests):
    index={i:i for i in range(len(bests))}
    fil1,fil2={},{}
    fils=[]
    taille=10*len(bests)
    while len(index)>1:
        for n in range(taille//2):
            i=random.choice(list(index.values()))
            index.pop(i)
            j=random.choice(list(index.values()))
            coupe=random.randint(0,len(bests)-1)
            counter=0
            for dinucleotide in bests[i].getTable():
                if counter<=coupe:
                    fil1[dinucleotide]=bests[i].getTable()[dinucleotide]
                    fil2[dinucleotide]=bests[j].getTable()[dinucleotide]
                if counter>coupe:
                    fil1[dinucleotide]=bests[j].getTable()[dinucleotide]
                    fil2[dinucleotide]=bests[i].getTable()[dinucleotide]
                counter+=1
            a=RotTable()
            a.newTable(fil1)
            fils.append(a)
            a.newTable(fil2)
            fils.append(a)
            fil1,fil2={},{}
            index[i]=i
        index.pop(i)
        index.pop(j)
    return fils



def crossover(bests):
    index1={i:i for i in range(len(bests)//2)}
    index2={j:j for j in range(len(bests))}
    fil={}
    fils=[]
    n_fils=20
    for k in range(len(bests)//2):
    #while len(index1)>0:
        i=random.choice(list(index1.values()))
        index2.pop(i)
        for n in range(n_fils):
            j=random.choice(list(index2.values()))
            for dinucleotide in bests[i].getTable():
                fil[dinucleotide]=random.choice([bests[i].getTable()[dinucleotide],bests[j].getTable()[dinucleotide]])
            a=RotTable()
            a.newTable(fil)
            fils.append(a)
            fil={}
        index2[i]=i
        index1.pop(i)
    return fils


def crossover2(bests):
    index={i:i for i in range(len(bests))}
    fil={}
    fils=[]
    taille=10*len(bests)
    while len(index)>1:
        for n in range(taille):
            i=random.choice(list(index.values()))
            index.pop(i)
            j=random.choice(list(index.values()))
            coupe=random.randint(0,len(bests)-1)
            counter=0
            for dinucleotide in bests[i].getTable():
                if counter<=coupe:
                    fil[dinucleotide]=bests[i].getTable()[dinucleotide]
                if counter>coupe:
                    fil[dinucleotide]=bests[j].getTable()[dinucleotide]
                counter+=1
            a=RotTable()
            a.newTable(fil)
            fils.append(a)
            fil={}
            index[i]=i
        index.pop(i)
        index.pop(j)
    return fils        


def crossover3(bests):
    index={i:i for i in range(len(bests))}
    fil1,fil2={},{}
    fils=[]
    taille=10*len(bests)
    while len(index)>1:
        for n in range(taille//2):
            i=random.choice(list(index.values()))
            index.pop(i)
            j=random.choice(list(index.values()))
            coupe=random.randint(0,len(bests)-1)
            counter=0
            for dinucleotide in bests[i].getTable():
                if counter<=coupe:
                    fil1[dinucleotide]=bests[i].getTable()[dinucleotide]
                    fil2[dinucleotide]=bests[j].getTable()[dinucleotide]
                if counter>coupe:
                    fil1[dinucleotide]=bests[j].getTable()[dinucleotide]
                    fil2[dinucleotide]=bests[i].getTable()[dinucleotide]
                counter+=1
            a=RotTable()
            a.newTable(fil1)
            fils.append(a)
            a.newTable(fil2)
            fils.append(a)
            fil1,fil2={},{}
            index[i]=i
        index.pop(i)
        index.pop(j)
    return fils


{'AA': [35.61792936393656, 7.922530158825857, -154], 'AC': [36.57921575262207, 6.629795068144642, 143], 'AG': [27.086778097291866, 5.980880915406881, 2], 'AT': [32.83175884418437, -0.8387458915322732, 0], 'CA': [34.27513350767127, -42.01809274689457, -64], 'CC': [33.75197175178809, 1.4782613664011999, -57], 'CG': [30.51738850516329, 7.106996926706735, 0], 'CT': [26.555009579857575, 4.811053093490184, -2], 'GA': [37.99202561701745, 7.684225796850326, 120], 'GC': [39.332809708627664, 4.173895658831757, 180], 'GG': [33.66079179641344, 4.228926463685514, 57], 'GT': [33.575351517535054, -4.948331239975031, -143], 'TA': [37.43055792309639, -0.43506698811919164, 0], 'TC': [35.29101019356463, -2.058702636563572, -120], 'TG': [33.56204789927513, -12.082552339606432, 64], 'TT': [35.643660139109734, 6.726585304530609, 154]}













class GA:
    """
    Classe permettant d'éxécuter un algorithme génétique sur un problème d'optimisation.
    """
    def genesis(n):
        #n is the number of individuals. Creates n RotTable objects
        pob=[RotTable() for i in range(n)]
        return pob

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

    def evaluerGA():
        pop = initialisation_population()
        while critere:  # Critere has to be defined
            va = pop.evaluation()

            pop1 = pop.selection(pop, va)
            pop1 = pop.mutation()
            pop = pop1

        return pop

