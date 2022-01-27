import sys  
from pathlib import Path  
import time
file = Path(__file__).resolve()  
package_root_directory = file.parents[1]  
sys.path.append(str(package_root_directory))

import matplotlib.pyplot as plt
import argparse
from os import popen
from re import A
from Algo_gen.Population import *

parser = argparse.ArgumentParser()
parser.add_argument("--filename", help="input filename of DNA sequence")
parser.parse_args()
args = parser.parse_args()

class Bench_mark:
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
        distance_vector=[]
        time_vector=[]
        iteration_vector=[]
        current_distance=math.inf       # Current_distance between the first and last elements of the chromosome 
        current_iteration=0                     # Iteration represents the current interaction
        population=Population(self.__pop_size)
        while current_distance>self.__max_dist and current_iteration<self.__max_iter:
            current_iteration+=1
            population.select_bests()                       # Selects the best 10% individuals in the population
            population.do_gen()                             # Makes a new população with crossovers between the best 10% individuals
            population.mutation()                           # Makes mutations in the population
            population.fitness(self.__seq)                  # Evalue the population
            best,current_distance=population.pickbest()     # Gets the best element of the population and the distance between its first and last gene 
            t=time.time()
            print("Iteration",current_iteration)
            print("Distance :",current_distance)    
            time_vector.append(t)
            distance_vector.append(current_distance)
            iteration_vector.append(current_iteration)

        return time_vector,distance_vector,iteration_vector

def plots():

    archive = open("seq_512.fasta",'r')
    for seq_512 in archive:
        seq_512 = seq_512.strip()
    archive.close()

    archive = open("plasmid_8k.fasta",'r')
    for seq_8k in archive:
        seq_8k = seq_8k.strip()
    archive.close()

    archive = open("plasmid_180k.fasta",'r')
    for seq_180k in archive:
        seq_180k = seq_180k.strip()
    archive.close()

    benchmark=Bench_mark(20,100,seq_512,0.01,args)
    time_512,distance_512,iteration_512=benchmark.darwin()

    benchmark=Bench_mark(20,100,seq_8k,0.01,args)
    time_8k,distance_8k,iteration_8k=benchmark.darwin()

    benchmark=Bench_mark(20,100,seq_180k,0.01,args)
    time_180k,distance_180k,iteration_180k=benchmark.darwin()

    fig = plt.figure()
    plt.xlabel('Iteration')
    plt.ylabel('Best distance')
    plt.title('Best distance x Iteration')
    plt.plot(iteration_512,distance_512,'b',iteration_8k,distance_8k,'r',iteration_180k,distance_180k,'g')
    plt.legend(['512 dinucleotides',"8k nucleotides","180k dinucleotides"])
    plt.savefig("Graph_distanceXiteration.png",format='png')
    plt.show()

    fig = plt.figure()
    plt.xlabel('Time(s)')
    plt.ylabel('Best distance')
    plt.title('Best distance x Time')
    plt.plot(time_512,distance_512,'b',time_8k,distance_8k,'r',time_180k,distance_180k,'g')
    plt.legend(['512 dinucleotides',"8k nucleotides","180k dinucleotides"])
    plt.savefig("Graph_distanceXtime.png",format='png')
    plt.show()






# Command pour running with a file : python main.py --filename=plasmid_8k.fasta

def main():
    # Running of genetics algorithm
    plots()

    # Running the Monte Carlo algorithm
    #else p==2: 
        #seq="AAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCCAGTAAACGAAAAAACCGCCTGGGGAGGCGGTTTAGTCGAAGGTTAAGTCAGTTGGGGACTGCTTAACCGGGTAACTGGCTTGGTGGAGCACAGATACCAAATACTGTCCTTCTAGTGTAGCCGCAGTTAGGCCACCACTTCAAGAACTCTTAATATCTCAATCCACCTTGTCCAGTTACCAGTGGCTGCTGCCAGTGGCGCTTTGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAGCCGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTATGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAACGGCGCAGCCTTTTCCTGGTTCTCGTTTTTTGCTCACATGTTTCTTTTGGCGTTATCCCCTGATTCTGTGGATAACCGCATCTCCGCTTTTGAGTGAGCAGACACCGCTCGCCGCAGCCGAACGACCGAGTGTAGCGAGTCAGTGAGCGAGGAAGCGGAAGAGCGCCGGAACGTGCATTTTCTCCTTACGCATCTGTGCGGCATTTCACATCGGACATGGTGCGCTTTCCATACAATTCGTACTGATGCCGCATAGTTAAGCCAGTATACACTCCGCTATCGCTACGTGACTGGTTCAGGGCTTCGCCCCGAAACCCCCTGACGCGCCCTGAGGGGCTTGTCTGCTCCCGGCATCCGCTCACAGACAAGCTGTTACCGTCTCCGGGAGCTGTATGTGTCAGAGGTTTTCACCGTCATCCCCGAAGCGTGCGA"
        #algo=GA(40,100,seq,5,args)
        #algo.darwin()


if __name__ == "__main__" :
    main()