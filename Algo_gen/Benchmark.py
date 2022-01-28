import sys  
from pathlib import Path  
import time
file = Path(__file__).resolve()  
package_root_directory = file.parents[1]  
sys.path.append(str(package_root_directory))

import matplotlib.pyplot as plt
from os import popen
from re import A
from Algo_gen.Population import *

class Benchmark:
    def __init__(self,size,max_iter,seq,dist_max):
        self.__pop_size = size          # size is the population size
        self.__max_iter=max_iter        # max_iter is the maximum number of interactions
        self.__max_dist=dist_max        # dist_max is the maximum distance acceptable    
        self.__seq=seq                  # seq is a gene sequence
    
    def darwin(self):
        t0=time.time()
        distance_vector=[]
        time_vector=[]
        iteration_vector=[]
        current_distance=math.inf               # Current_distance between the first and last elements of the chromosome 
        current_iteration=0                     # Iteration represents the current interaction
        population=Population(self.__pop_size)
        while current_distance>self.__max_dist and current_iteration<self.__max_iter:
            current_iteration+=1
            population.select_bests()                       # Selects the best 10% individuals in the population
            population.do_gen()                             # Makes a new população with crossovers between the best 10% individuals
            population.mutation()                           # Makes mutations in the population
            population.fitness(self.__seq)                  # Evalue the population
            best,current_distance=population.pickbest()     # Gets the best element of the population and the distance between its first and last gene 
            t1=time.time()
            print("Iteration",current_iteration)
            print("Distance :",current_distance)    
            time_vector.append(t1-t0)
            distance_vector.append(current_distance)
            iteration_vector.append(current_iteration)


        traj=Traj3D()
        traj.compute(self.__seq,best)
        # Save of each solution
        if len(self.__seq)<2000:
            archive = open("seq_512_solution_AG.fasta",'w')     
            archive.write(str(best.get_chr()))
            archive.close()
            traj.compute(self.__seq,best)
            traj.draw("seq_512_solution_AG.png")

        elif len(self.__seq)>2000 and len(self.__seq)<10000:
            archive = open("seq_8k_solution_AG.fasta",'w')     
            archive.write(str(best.get_chr()))
            archive.close()
            traj.compute(self.__seq,best)
            traj.draw("seq_8k_solution_AG.png")
        
        elif len(self.__seq)>10000:
            archive = open("seq_180k_solution_AG.fasta",'w')   
            archive.write(str(best.get_chr()))
            archive.close()
            traj.compute(self.__seq,best)
            traj.draw("seq_180k_solution_AG.png")


        return time_vector,distance_vector,iteration_vector

def plot_benchmark():
    # Reading the sequences
    archive = open("plasmid_512.fasta",'r')
    for seq_512 in archive:
        seq_512 = seq_512.strip()
    archive.close()

    archive = open("plasmid_8k.fasta",'r')
    lineList = [line.rstrip('\n') for line in archive]
    seq_8k = ''.join(lineList[1:])
    archive.close()

    archive = open("plasmid_180k.fasta",'r')
    lineList = [line.rstrip('\n') for line in archive]
    seq_180k = ''.join(lineList[1:])
    archive.close()

    # It makes a benchmark test for each sequences
    benchmark_512=Benchmark(40,250,seq_512,0.01)
    time_512,distance_512,iteration_512=benchmark_512.darwin()

    benchmark_8k=Benchmark(40,250,seq_8k,0.01)
    time_8k,distance_8k,iteration_8k=benchmark_8k.darwin()

    benchmark_180k=Benchmark(40,250,seq_180k,0.01)
    time_180k,distance_180k,iteration_180k=benchmark_180k.darwin()

    # Graph_distance X Iteration plot
    fig = plt.figure()
    plt.xlabel('Iteration')
    plt.ylabel('Best distance')
    plt.title('Best distance x Iteration')
    plt.plot(iteration_512,distance_512,'b',iteration_8k,distance_8k,'r',iteration_180k,distance_180k,'g')
    plt.legend(['512 dinucleotides',"8k nucleotides","180k dinucleotides"])
    plt.savefig("Graph_distanceXiteration_AG.png",format='png')
    plt.show()

    # Graph_distance X time plot
    fig = plt.figure()
    plt.xlabel('Time(s)')
    plt.ylabel('Best distance')
    plt.title('Best distance x Time')
    plt.plot(time_512,distance_512,'b',time_8k,distance_8k,'r',time_180k,distance_180k,'g')
    plt.legend(['512 dinucleotides',"8k dinucleotides","180k dinucleotides"])
    plt.savefig("Graph_distanceXtime_AG.png",format='png')
    plt.show()



def main():
    plot_benchmark()


if __name__ == "__main__" :
    main()