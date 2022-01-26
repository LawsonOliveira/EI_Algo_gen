
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