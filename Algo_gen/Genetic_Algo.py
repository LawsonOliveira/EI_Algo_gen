import sys  
from pathlib import Path  
file = Path(__file__).resolve()  
package_root_directory = file.parents[1]  
sys.path.append(str(package_root_directory))

from os import popen
from re import A
from Algo_gen.Population import *

# This class represents genetics algorithm
class Genetic_Algo:
    def __init__(self,size,max_iter,seq,dist_max,args=None):
                                        # seq is a gene sequence
        self.__pop_size = size          # size is the population size
        self.__max_iter=max_iter        # max_iter is the maximum number of interactions
        self.__args=args                # args is the argparse arguments
        self.__max_dist=dist_max        # dist_max is the maximum distance acceptable    
        if self.__args.filename:
            # Read file
            lineList = [line.rstrip('\n') for line in open(self.__args.filename)]
            # Formatting
            self.__seq = ''.join(lineList[1:])
        else:
            self.__seq=seq


    def darwin(self):
        current_distance=math.inf       # Current_distance between the first and last elements of the sequence 
        iteration=0                     # Iteration represents the current interaction
        population=Population(self.__pop_size)
        while current_distance>self.__max_dist and iteration<self.__max_iter:
            iteration+=1
            population.select_bests()                       # Selects the best 10% individuals in the population
            population.do_gen()                             # Makes a new população with crossovers between the best 10% individuals
            population.mutation()                           # Makes mutations in the population
            population.fitness(self.__seq)                  # Evalue the population
            best,current_distance=population.pickbest()     # Gets the best element of the population and the distance between its first and last gene 
            print("Iteration",iteration)
            print("Distance :",current_distance)


        traj=Traj3D()
        if self.__args.filename:                            # Graph plot of a input file in the terminal
            archive = open(self.__args.filename+"_solution_.fasta",'w')
            archive.write(str(best.get_chr()))
            archive.close()
            traj.compute(self.__seq,best)
            traj.draw(self.__args.filename+".png")

        else:                                               # Graph plot of a sequence in the main fonction
            archive = open("exemple"+"_solution_.fasta",'w')
            archive.write(str(best.get_chr()))
            archive.close()
            traj.compute(self.__seq,best)
            traj.draw("sample.png")

