import sys  
from pathlib import Path  
file = Path(__file__).resolve()  
package_root_directory = file.parents[1]  
sys.path.append(str(package_root_directory))

from Algo_gen.Chromosome import *
import numpy as np
from Traj3D import *
import random

class Population:
    
    def __init__(self,n,size_dinucle=16):
        self.__pop_size=n
        self.__distanceinucle_size=size_dinucle
        self.__pop = [Chromosome() for i in range(n)]
        self.__seq=[]
        self.__distance=[math.inf for i in range(n)]
        self.__bests=[None for i in range(n//10)]

    def __str__(self):
        s=[]
        for i in self.__pop:
            s.append(str(i))
        return str(s)

    def fitness(self,seq):
        self.__seq=seq
        self.__distance=np.array([math.inf for k in range(self.__pop_size)]) 
        trajs=[Traj3D() for i in range(self.__pop_size)]         #initialize the trajs                    #initialize the vector that has the distances
        for j in range(self.__pop_size):
            trajs[j].compute(self.__seq,self.__pop[j])        #calculates each traj
            xyz = np.array(trajs[j].getTraj())  
            x, y, z = xyz[:,0], xyz[:,1], xyz[:,2]
            self.__distance[j]=np.sqrt((x[-1])**2+(y[-1])**2+(z[-1])**2)+0.01*(np.dot(xyz[1],xyz[-1]+xyz[-2])+np.linalg.norm(xyz[1])*np.linalg.norm(xyz[-1]+xyz[-2]))#calculates the distance from the tip of the chain to the center

    def select_bests(self):
        #Chooses the better 10% of the population
        for j in range(self.__pop_size//10):
            i=np.argmin(self.__distance)
            self.__distance[i]=math.inf
            self.__bests[j]=self.__pop[i]


    def do_gen(self):
        index = [[random.randint(0,len(self.__bests)-1),random.randint(0,len(self.__bests)-1)] for i in range(int(len(self.__bests)*10//2))]
        new_gen=[]
        for i,j in index:
            gene=Chromosome()
            enfant1,enfant2=gene.crossover(self.__bests[i],self.__bests[j])
            new_gen.append(Chromosome(enfant1))
            new_gen.append(Chromosome(enfant2))
        self.__pop=new_gen
        return self
            
    def pickbest(self):
        #Picks the best individual from pob
        i=np.argmin(self.__distance)
        return self.__pop[i],self.__distance[i]
    
    def mutation(self,taux=0.05):
        # Makes the mutation in the population
        numb_of_dinucleotides=self.__distanceinucle_size 
        n_mutates=int(self.__pop_size*numb_of_dinucleotides*taux)       # This is the total number of mutations there will be
        for i in range(n_mutates):
            individu=random.randint(0,len(self.__pop)-1)
            self.__pop[individu].apply_mutation()

    def new_pop(self,new_pop):
        # Update the population
        self.__pop=new_pop

    def add_individual(self, individual,i):
        # Adds a new element to the population at position i
        self.__pop[i]=individual
    
    def get_size(self):
        # Return the population size
        return self.__pop_size

    def get_indiv(self, i):
        # Return the individual in the position i of the population
        return self.__pop[i]

    def get_pop(self):
        # Return the population
        return self.__pop
    
    def get_distance(self):
        # Return the distance between the first gene and last gene of each chromosome
        return self.__distance

    def get_bests(self):
        # Return the bests individuals of the population
        return self.__bests