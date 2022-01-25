
import random

table1 = {
    "AA": [35.62, 7.2, -154,      0.06,  0.6, 0],
    "AC": [34.4, 1.1,  143,      1.3,  5, 0],
    "AG": [27.7, 8.4,    2,      1.5,  3, 0],
    "AT": [31.5, 2.6,    0,      1.1,  2, 0],
    "CA": [34.5, 3.5,  -64,      0.9, 34, 0],
    "CC": [33.67, 2.1,  -57,      0.07,  2.1, 0],
    "CG": [29.8, 6.7,    0,      1.1,  1.5, 0],
    "CT": [27.7, 8.4,   -2,      1.5,  3, 0],
    "GA": [36.9, 5.3,  120,      0.9,  6, 0],
    "GC": [40, 5,  180,      1.2,  1.275, 0],
    "GG": [33.67, 2.1,   57,      0.07,  2.1, 0],
    "GT": [34.4, 1.1, -143,      1.3,  5, 0],
    "TA": [36, 0.9,    0,      1.1,  2, 0],
    "TC": [36.9, 5.3, -120,      0.9,  6, 0],
    "TG": [34.5, 3.5,   64,      0.9, 34, 0],
    "TT": [35.62, 7.2,  154,      0.06,  0.6, 0]
}
table2 = {
    "AA": [50, 7.2, -154,      0.06,  0.6, 0],
    "AC": [34.4, 1.1,  143,      1.3,  5, 0],
    "AG": [50.7, 8.4,    2,      1.5,  3, 0],
    "AT": [31.5, 2.6,    0,      1.1,  2, 0],
    "CA": [50.5, 3.5,  -64,      0.9, 34, 0],
    "CC": [33.67, 2.1,  -57,      0.07,  2.1, 0],
    "CG": [50.8, 6.7,    0,      1.1,  1.5, 0],
    "CT": [27.7, 8.4,   -2,      1.5,  3, 0],
    "GA": [36.9, 5.3,  120,      0.9,  6, 0],
    "GC": [50, 5,  180,      1.2,  1.275, 0],
    "GG": [33.67, 2.1,   57,      0.07,  2.1, 0],
    "GT": [50.4, 1.1, -143,      1.3,  5, 0],
    "TA": [36, 0.9,    0,      1.1,  2, 0],
    "TC": [50.9, 5.3, -120,      0.9,  6, 0],
    "TG": [34.5, 3.5,   64,      0.9, 34, 0],
    "TT": [50.62, 7.2,  154,      0.06,  0.6, 0]
}

all_nuclets=[table1,table2,table1,table2]

def crossover(bests):
    index={i:i for i in range(len(all_nuclets))}
    fil={}
    fils=[]
    taille=10
    while len(index)>0:
        i=random.choice(index)
        index.pop(i)
        j=random.choice(index)
        index.pop(j)
        print(index)
        for n in range(taille):
            for dinucleotide_pere in all_nuclets[i]:
                for dinucleotide_mere in all_nuclets[j]:
                    fil[dinucleotide_pere]=random.choice([all_nuclets[i][dinucleotide_pere],all_nuclets[j][dinucleotide_mere]])
            fils.append(newTable(fil))
            fil={}
    return fils
        
print(crossover(all_nuclets))
