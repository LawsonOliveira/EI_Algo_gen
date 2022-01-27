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
    print(random.random)
    if random.random() < p:
        return 1
    else:
        return 0

def genesis(n):
    #n is the number of individuals. Creates n RotTable objects
    pob=[RotTable() for i in range(n)]
    return pob
                
def mutation(pob):
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


def crossover2(bests):
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
        aux=RotTable()
        aux.newTable(fil1)
        fils.append(aux)
        aux.newTable(fil2)
        fils.append(aux)
        fil1,fil2={},{}
    return fils

def crossover1(pere,mere,prob,keys):
    enfant1,enfant2={},{}
    pos=0
    for key in keys:
        if prob[pos]<0.5:
            enfant1[key]=pere[pos]
            enfant2[key]=mere[pos]
        else:
            enfant1[key]=mere[pos]
            enfant2[key]=pere[pos]
        pos+=1
    return enfant1,enfant2

def crossover(pere,mere):
    c=random.randint(0,len(pere)-1)
    enfant1=np.append(pere[:c],mere[c:])
    enfant2=np.append(pere[c:],mere[:c])
    return enfant1,enfant2

def reproduction(bests):
    index = [[random.randint(0,len(bests)-1),random.randint(0,len(bests)-1)] for i in range(int(len(bests)*10//2))]
    enfants=[]
    for i,j in index:
        enfant1,enfant2=crossover(copy.deepcopy(bests[i].getTable),bests[j])
        table=RotTable()
        table.newTable(enfant1)
        enfants.append(table)
        table.newTable(enfant2)
        enfants.append(table)
    return enfants



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
        pob=reproduction(pob)
        t1=time.time()
        print(t1-t0)
        pob=mutation(pob)
        t2=time.time()

        print(t2-t1)
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
        best={'AA': [35.61792936393656, 7.922530158825857, -154], 'AC': [36.57921575262207, 6.629795068144642, 143], 'AG': [27.086778097291866, 5.980880915406881, 2], 'AT': [32.83175884418437, -0.8387458915322732, 0], 'CA': [34.27513350767127, -42.01809274689457, -64], 'CC': [33.75197175178809, 1.4782613664011999, -57], 'CG': [30.51738850516329, 7.106996926706735, 0], 'CT': [26.555009579857575, 4.811053093490184, -2], 'GA': [37.99202561701745, 7.684225796850326, 120], 'GC': [39.332809708627664, 4.173895658831757, 180], 'GG': [33.66079179641344, 4.228926463685514, 57], 'GT': [33.575351517535054, -4.948331239975031, -143], 'TA': [37.43055792309639, -0.43506698811919164, 0], 'TC': [35.29101019356463, -2.058702636563572, -120], 'TG': [33.56204789927513, -12.082552339606432, 64], 'TT': [35.643660139109734, 6.726585304530609, 154]}
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
