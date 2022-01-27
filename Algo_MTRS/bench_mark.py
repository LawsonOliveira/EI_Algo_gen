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
from Algo_MTRS.montecarlo_fonctionarbre import *


class Bench_mark:
    def __init__(self,nbit,seq):
                                        # seq is a gene sequence
        self.__nbit = nbit              # nbit       
        self.__seq=seq


    def compute(self,critere=10**(-3)):
        root=node()
        time_vector=[]
        distance_vector=[]
        iteration_vector=[]
        t_start=time.time()
        root.writeValeur(-np.inf)
        value = critere+1  # Value of the first root, MUST NOT BE 0 WATCH OUT
        nbiteration = 0
        print("initialisation", value)
        print()
        while nbiteration < self.__nbit and abs(value) > critere:
            backpropagation(root,self.__seq)  # We backpropagation like always
            print(" nbiteration", nbiteration)
            print("value", value)
            
            nbiteration += 1
            value = root.getvalue()

            distance_vector.append(-value)
            time_vector.append(time.time()-t_start)
            iteration_vector.append(nbiteration)
        

        return time_vector,distance_vector,iteration_vector

def plots():

    archive = open("seq_512.fasta",'r')
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

    benchmark=Bench_mark(500,seq_512)
    time_512,distance_512,iteration_512=benchmark.compute()

    benchmark=Bench_mark(500,seq_8k)
    time_8k,distance_8k,iteration_8k=benchmark.compute()

    benchmark=Bench_mark(500,seq_180k)
    time_180k,distance_180k,iteration_180k=benchmark.compute()

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

    
if __name__ == "__main__" :
    main()