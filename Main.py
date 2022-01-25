from RotTable import *
from Traj3D import *
import random

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--filename", help="input filename of DNA sequence")
parser.parse_args()
args = parser.parse_args()

def bernoulli(p):
    if random.random() < p:
        return 1
    else:
        return 0

def genesis(n):
    #n is the number of individuals. Creates n RotTable objects
    pob=[RotTable() for i in range(n)]
    return pob


def mutation(pob):
    #Applies mutations tu the population
    p=0.05 #change this equation please
    for i in pob:
        if bernoulli(p)==1:
            for j in range(3): #three is the number of mutations. Can be changed
                u=random.randint(0,15)
                i.mut(u)       #Applies a mutation to a single individual in a single nucleotid
    return pob
                
def selection(pob,D,n):
    #Chooses the better half of the population
    new_pob=[]
    for j in range(int(n/10)):
        i=np.argmin(D)
        new_pob.append(pob[i])
    return new_pob

def crossover(bests):
    index={i:i for i in range(len(bests))}
    fil={}
    fils=[]
    taille=10*len(bests)
    while len(index)>1:
        i=random.choice(list(index.values()))
        index.pop(i)
        j=random.choice(list(index.values()))
        index.pop(j)
        for n in range(taille):
            for dinucleotide_pere in bests[i].getTable():
                for dinucleotide_mere in bests[j].getTable():
                    fil[dinucleotide_pere]=random.choice([bests[i].getTable()[dinucleotide_pere],bests[j].getTable()[dinucleotide_mere]])
            a=RotTable()
            a.writeTable(fil)
            fils.append(a)
            fil={}
    return fils
        
def pickbest(pob,D):
    #Picks the best individual from pob
    i=np.argmin(D)
    return pob[i]

def darwin(pob,n,k,seq):
    #takes a poblation and a the number of repetitions k
    trajs=[Traj3D() for i in range(n)]         #initialize the trajs
    D=[0 for m in range(n)]                    #initialize the vector that has the distances

    for i in range(k):
        for j in range(n):
            trajs[j].compute(seq,pob[j])        #calculates each traj
            xyz = np.array(trajs[j].getTraj())  
            x, y, z = xyz[:,0], xyz[:,1], xyz[:,2]
            D[j]=np.sqrt(x[-1]**2 + y[-1]**2 + z[-1]**2 ) #calculates the distance from the tip of the chain to the center
            
        pob=selection(pob,D,n)
        pob=crossover(pob)
        pob=mutation(pob)




    best=pickbest(pob,D)
    return best
        

def main():
    n=100
    k=10
    pob = genesis(n)
    traj = Traj3D()

    if args.filename:
        # Read file
        lineList = [line.rstrip('\n') for line in open(args.filename)]
        # Formatting
        seq = ''.join(lineList[1:])
        best_ind=darwin(pob,n,k,seq)
        traj.compute(seq,best_ind)
    else:
        best_ind=darwin(pob,n,k,"AAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCCAGTAAACGAAAAAACCGCCTGGGGAGGCGGTTTAGTCGAAGGTTAAGTCAGAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCCAGTAAACGAAAAAACCGCCTGGGGAGGCGGTTTAGTCGAAGGTTAAGTCAG")
        traj.compute("AAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCCAGTAAACGAAAAAACCGCCTGGGGAGGCGGTTTAGTCGAAGGTTAAGTCAGAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCCAGTAAACGAAAAAACCGCCTGGGGAGGCGGTTTAGTCGAAGGTTAAGTCAG", best_ind)

    print(traj.getTraj())

    if args.filename:
        traj.draw(args.filename+".png")
    else:
        traj.draw("sample.png")


if __name__ == "__main__" :
    main()
