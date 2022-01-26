from RotTable import *
from Traj3D import *
import random
import argparse
import time
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
    p=0.05 #change this equation please
    for i in pob:
        for k in range(len(i.getTable())):
            if bernoulli(p)==1:
                i.mut(k)       #Applies a mutation to a single individual in a single nucleotid
    return pob

def mutation2(pob):
    numb_of_dinucleotides=len(pob[0].getTable())
    p=0.05 #change this equation please
    n_mutates=int(len(pob)*numb_of_dinucleotides*p)
    for i in range(n_mutates):
        individu=random.randint(0,len(pob)-1)
        dinucleotide=random.randint(0,numb_of_dinucleotides-1)
        pob[individu].mut(dinucleotide)
    return pob


def selection(pob,D,n):
    #Chooses the better 10% of the population
    new_pob=[]
    for j in range(n//10):
        i=np.argmin(D)
        D[i]=math.inf
        new_pob.append(pob[i])
    return new_pob

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

def crossover4(bests):
    index1=[int(2*i) for i in range(len(bests)//2)]
    index2=[int(2*i+1) for i in range(len(bests)//2)]
    fil1,fil2={},{}
    fils=[]
    taille=int(10*len(bests))
    for n in range(taille//2):
        i=random.choice(index1)
        j=random.choice(index2)
        counter=0
        coupe=random.randint(0,len(bests)-1)
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
    return fils


def pickbest(pob,D):
    #Picks the best individual from pob
    i=np.argmin(D)
    return pob[i],D[i]

def fitness(pob,n,D,seq):
    trajs=[Traj3D() for i in range(n)]         #initialize the trajs                    #initialize the vector that has the distances
    for j in range(n):
        trajs[j].compute(seq,pob[j])        #calculates each traj
        xyz = np.array(trajs[j].getTraj())  
        x, y, z = xyz[:,0], xyz[:,1], xyz[:,2]
        #D[j]=np.sqrt((x[-1])**2+(y[-1])**2+(z[-1])**2)#calculates the distance from the tip of the chain to the center
        D[j]=np.sqrt((x[-1])**2+(y[-1])**2+(z[-1])**2)+0.01*(np.dot(xyz[1],xyz[-1]+xyz[-2])+np.linalg.norm(xyz[1])*np.linalg.norm(xyz[-1]+xyz[-2]))#calculates the distance from the tip of the chain to the center

    return D
    


def darwin(pob,n,k,seq):
    #takes a poblation and a the number of repetitions k
    D=[100 for k in range(n)]        
    p=0
    distance_actuel=100
    while distance_actuel>1 and p<k:
        #for i in range(k):
        pob=selection(pob,D,n)
        t0=time.time()
        pob=crossover4(pob)
        t1=time.time()
        print(t1-t0)
        pob=mutation2(pob)
        D=fitness(pob,n,copy.deepcopy(D),seq)
        best,distance_actuel=pickbest(pob,D)
        p+=1
        print("Iteration",p)
        print("Distance minimale",distance_actuel)
    return best
  
def main():
    n=400                #n>40 et n//4=0; k>200
    k=100
    pob = genesis(n)
    traj=Traj3D()
    if args.filename:
        # Read file
        lineList = [line.rstrip('\n') for line in open(args.filename)]
        # Formatting
        seq = ''.join(lineList[1:])
        best=darwin(pob,n,k,seq)
        archive = open(args.filename+"_solution_.fasta",'w')
        archive.write(str(best.getTable()))
        archive.close()
        traj.compute(seq,best)
        traj.draw(args.filename+".png")
    else:
        best=darwin(pob,n,k,"AAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCCAGTAAACGAAAAAACCGCCTGGGGAGGCGGTTTAGTCGAAGGTTAAGTCAGTTGGGGACTGCTTAACCGGGTAACTGGCTTGGTGGAGCACAGATACCAAATACTGTCCTTCTAGTGTAGCCGCAGTTAGGCCACCACTTCAAGAACTCTTAATATCTCAATCCACCTTGTCCAGTTACCAGTGGCTGCTGCCAGTGGCGCTTTGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAGCCGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTATGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAACGGCGCAGCCTTTTCCTGGTTCTCGTTTTTTGCTCACATGTTTCTTTTGGCGTTATCCCCTGATTCTGTGGATAACCGCATCTCCGCTTTTGAGTGAGCAGACACCGCTCGCCGCAGCCGAACGACCGAGTGTAGCGAGTCAGTGAGCGAGGAAGCGGAAGAGCGCCGGAACGTGCATTTTCTCCTTACGCATCTGTGCGGCATTTCACATCGGACATGGTGCGCTTTCCATACAATTCGTACTGATGCCGCATAGTTAAGCCAGTATACACTCCGCTATCGCTACGTGACTGGTTCAGGGCTTCGCCCCGAAACCCCCTGACGCGCCCTGAGGGGCTTGTCTGCTCCCGGCATCCGCTCACAGACAAGCTGTTACCGTCTCCGGGAGCTGTATGTGTCAGAGGTTTTCACCGTCATCCCCGAAGCGTGCGA")
        archive = open("exemple"+"_solution_.fasta",'w')
        archive.write(str(best.getTable()))
        archive.close()
        traj.compute("AAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCCAGTAAACGAAAAAACCGCCTGGGGAGGCGGTTTAGTCGAAGGTTAAGTCAGTTGGGGACTGCTTAACCGGGTAACTGGCTTGGTGGAGCACAGATACCAAATACTGTCCTTCTAGTGTAGCCGCAGTTAGGCCACCACTTCAAGAACTCTTAATATCTCAATCCACCTTGTCCAGTTACCAGTGGCTGCTGCCAGTGGCGCTTTGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAGCCGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTATGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAACGGCGCAGCCTTTTCCTGGTTCTCGTTTTTTGCTCACATGTTTCTTTTGGCGTTATCCCCTGATTCTGTGGATAACCGCATCTCCGCTTTTGAGTGAGCAGACACCGCTCGCCGCAGCCGAACGACCGAGTGTAGCGAGTCAGTGAGCGAGGAAGCGGAAGAGCGCCGGAACGTGCATTTTCTCCTTACGCATCTGTGCGGCATTTCACATCGGACATGGTGCGCTTTCCATACAATTCGTACTGATGCCGCATAGTTAAGCCAGTATACACTCCGCTATCGCTACGTGACTGGTTCAGGGCTTCGCCCCGAAACCCCCTGACGCGCCCTGAGGGGCTTGTCTGCTCCCGGCATCCGCTCACAGACAAGCTGTTACCGTCTCCGGGAGCTGTATGTGTCAGAGGTTTTCACCGTCATCCCCGAAGCGTGCGA",best)
        traj.draw("sample.png")

if __name__ == "__main__" :
    main()
