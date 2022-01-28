import sys  
from pathlib import Path  
import time
file = Path(__file__).resolve()  
package_root_directory = file.parents[1]  
sys.path.append(str(package_root_directory))

from Algo_MCTS.MCTS import MCTS
import matplotlib.pyplot as plt
from os import popen
from re import A
from Algo_gen.Population import *
from Algo_MCTS.Monte_Carlo_fonctions import *


class Benchmark_MCTS:
    def __init__(self,nbit,seq):
                                        # seq is a gene sequence
        self.__nbit = nbit              # nbit       
        self.__seq=seq


    def compute(self,critere=10**(-3)):
        root=Node()
        time_vector=[]
        distance_vector=[]
        iteration_vector=[]
        t_start=time.time()
        root.writeValeur(-np.inf)
        value = critere+1  # Value of the first root, MUST NOT BE 0 WATCH OUT
        nbiteration = 0
        while nbiteration < self.__nbit and abs(value) > critere:
            backpropagation(root,self.__seq)  # We backpropagation like always
            print("Iteration:", nbiteration)
            print("Distance:", value)
            
            nbiteration += 1
            value = root.getvalue()

            distance_vector.append(-value)
            time_vector.append(time.time()-t_start)
            iteration_vector.append(nbiteration)
        return time_vector,distance_vector,iteration_vector

def plots():
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
    benchmark=Benchmark_MCTS(2500,seq_512)
    time_512,distance_512,iteration_512=benchmark.compute()

    benchmark=Benchmark_MCTS(1000,seq_8k)
    time_8k,distance_8k,iteration_8k=benchmark.compute()

    benchmark=Benchmark_MCTS(500,seq_180k)
    time_180k,distance_180k,iteration_180k=benchmark.compute()

    # Graph_distance X Iteration plot
    fig = plt.figure()
    plt.xlabel('Iteration')
    plt.ylabel('Best distance')
    plt.title('Best distance x Iteration')
    plt.plot(iteration_512,distance_512,'b',iteration_8k,distance_8k,'r',iteration_180k,distance_180k,'g')
    plt.legend(['512 nucleotides',"8k nucleotides","180k nucleotides"])
    plt.savefig("Graph_distanceXiteration_MCTS.png",format='png')
    plt.show()

    fig = plt.figure()
    plt.xlabel('Time(s)')
    plt.ylabel('Best distance')
    plt.title('Best distance x Time')
    plt.plot(time_512,distance_512,'b',time_8k,distance_8k,'r',time_180k,distance_180k,'g')
    plt.legend(['512 dinucleotides',"8k dinucleotides","180k dinucleotides"])
    plt.savefig("Graph_distanceXtime_MCTS.png",format='png')
    plt.show()

def main():
    plots()
    
if __name__ == "__main__" :
    main()
