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
    #n is the number of individuals
    pob=[RotTable() for i in range(n)]
    return pob


def mutation(pob):
    p=0.05 #change this equation please
    for i in pob:
        if bernoulli(p)==1:
            for j in range(3): #three is the number of mutations. Can be changed
                u=random.randint(0,15)
                i.mut(u)
    return pob
                
def selection(pob,D,n):
    new_pob=[]
    for j in range(int(n/2)):
        i=np.argmin(D)
        new_pob.append(pob[i])
    return new_pob

def crossover(pob):
    return pob

def pickbest(pob,D):
    i=np.argmin(D)
    return pob[i]

def darwin(pob,n,k,seq):
    #takes a poblation and a the number of repetitions k
    trajs=[Traj3D() for i in range(n)]
    for i in range(k):
        D=[0 for m in range(n)]
        for j in range(n):
            trajs[j].compute(seq,pob[j])
            xyz = np.array(trajs[j].getTraj())
            x, y, z = xyz[:,0], xyz[:,1], xyz[:,2]
            D[j]=np.sqrt(x[-1]**2 + y[-1]**2 + z[-1]**2 )
            
        pob=selection(pob,D,n)
        pob=crossover(pob)
        pob=mutation(pob)




    best=pickbest(pob,D)
    return best
        

def main():
    n=10
    k=1
    pob = genesis(n)
    traj = Traj3D()

    if args.filename:
        # Read file
        lineList = [line.rstrip('\n') for line in open(args.filename)]
        # Formatting
        seq = ''.join(lineList[1:])
        

        traj.compute(seq, rot_table)
    else:
        best_ind=darwin(pob,n,k,"AAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCCAGTAAACGAAAAAACCGCCTGGGGAGGCGGTTTAGTCGAAGGTTAAGTCAG")
        traj.compute("AAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCCAGTAAACGAAAAAACCGCCTGGGGAGGCGGTTTAGTCGAAGGTTAAGTCAG", best_ind)

    print(traj.getTraj())

    if args.filename:
        traj.draw(args.filename+".png")
    else:
        traj.draw("sample.png")


if __name__ == "__main__" :
    main()
