
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
