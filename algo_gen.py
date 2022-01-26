from os import popen
from re import A
from random import *
git 

class Gene():
    def __init__(self, bit, size):
        self.bit = bit
        self.size = size
    @classmethod
    def Generation(cls, size):
        #size: number of bits used to store the gene value
        bina = randint(0, 2**size-1)
        return Gene(bina, size)

    def Mutation(self,fct_mutation):
        fct_mutation(self)

    @classmethod
    def Crossover(cls, gene1, gene2, fct_crossover):
        return(fct_crossover(gene1, gene2))

class Chromosome():
    def __init__(self, genes=[], mutation_genes=[]):
        self.genes = genes
        self.mutation_genes = mutation_genes
    
    def nbr_genes(self):
        return len(self.genes)

    def add_gene(self, gene):
        self.genes.append(gene)
    
    def add_mutation_gene(self, mutation_gene):
        self.mutation_genes.append(mutation_gene)

    def __getitem__(self, i):
        return self.genes[i]

    def Mutation(self, moy):
        #moy: average number of bits to be muted in each gene
        for gene in self.genes:
            for i in range(gene.size):
                if random() < moy / gene.size:
                    gene.bin = 1 - gene.bin

    @classmethod
    def Crossover(cls, chromosome1, chromosome2, fct_crossover):
        assert chromosome1.nb_genes() == chromosome2.nb_genes()
        #we initialize chr_os1 and chr_os2: chromosome offsprings
        chr_os1 = Chromosome([], [])
        chr_os2 = Chromosome([], [])
        #crossover for each gene
        for i in range(chromosome1.nbr_genes()):
            g1, g2 = Gene.Croisement(
                chromosome1[i], chromosome2[i], fct_crossover)

            if chromosome1.mutation_genes != [] and chromosome2.mutation_genes != [] and len(chromosome1.genes_mutation) == len(chromosome2.genes_mutation):
                mu_gene1, mu_gene2 = Gene.Crossover(
                    chromosome1.mutation_genes[i], chromosome2.mutation_genes[i], fct_crossover)
                #Adding the new mutation genes to the offspring chromosomes
                chr_os1.add_mutation_gene(mu_gene1)
                chr_os2.add_mutation_gene(mu_gene2)
            
            chr_os1.add_gene(g1)
            chr_os2.add_gene(g2)

        return (chr_os1, chr_os2)



class Individual():
    def __init__(self, chromosomes=[]):
        self.chromosomes = chromosomes
        self.nbr_chr = len(chromosomes)

    def add_chromosome(self, chromosome):
        self.chromosomes.append(chromosome)

    def nbr_chromosomes(self):
        return len(self.chromosomes)

    def __getitem__(self, i):
        if type(i) == tuple:
            return(self.chromosomes[i[0]][i[1]])
        elif type(i) == int:
            return(self.chromosomes[i])

    def Mutation(self, moy):
        #moy: average number of bits to be muted in each gene
        for chromosome in self.chromosomes:
            chromosome.Mutation(moy)

    @classmethod
    def Crossover(cls, individu_1, individu_2, fct_crossover):
        #Initializing offspring individuals
        indiv_os1 = Individual([])
        indiv_os2 = Individual([])

        #Crossover for each chromosome
        for i in range(individual_1.nbr_chr):
            chr1, chr2 = Chromosome.Croisement(
                individual_1[i], individual_2[i], fct_crossover)
        
            indiv_os1.add_chromosome(c1)
            indiv_os2.add_chromosome(c2)

        return (indiv_os1, indiv_os2)




class population():

    #def __init__(self, selection, mutation, fit):
        #self.pop = []
        #self.selection = selection
        #self.mutation = mutation
        #self.fit = fit

    def __init__(self, individuals=[], rate=0.05):
        self.individuals = individuals
        self.taille = len(individuals)
        self.initial_size = len(individuals)
        self.size_best = round(len(individuals)*rate)

    def add_indiv(self, indiv):
        self.individus.append(indiv)
    
    def size(self):
        return len(self.individuals)

    def __getindiv__(self, i):
        return self.individuals[i]

    def Mutation(self, moy, fct_mutation):
        for indiv in self.individuals[self.size_best:]:
            Individual.Mutation(fct_mutation, moy)

    def Crossover(self, fct_crossover):
        new_individuals = []
        while (self.size() + len(new_individuals)) < self.initial_size:
            i, j = randint(0, self.size()-1), randint(0, self.size()-1)
            new_indiv1, new_indiv2 = Individual.Crossover(
                self.individuals[i], self.individuals[j], fct_crossover)
            new_individuals.extend([new_indiv1, new_indiv2])
        self.individuals += new_individuals
        if self.size() > self.initial_size:
            self.individuals.pop()



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

    def evaluerGA():
        pop = initialisation_population()
        while critere:  # Critere has to be defined
            va = pop.evaluation()

            pop1 = pop.selection(pop, va)
            pop1 = pop.mutation()
            pop = pop1

        return pop

