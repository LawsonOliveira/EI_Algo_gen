from Algo_gen.Chromosome import *
from Algo_gen.Genetic_Algo import *
from Algo_gen.Population import *
import numpy

population=Population(20)
seq="AAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCCAGTAAACGAAAAAACCGCCTGGGGAGGCGGTTTAGTCGAAGGTTAAGTCAGTTGGGGACTGCTTAACCGGGTAACTGGCTTGGTGGAGCACAGATACCAAATACTGTCCTTCTAGTGTAGCCGCAGTTAGGCCACCACTTCAAGAACTCTTAATATCTCAATCCACCTTGTCCAGTTACCAGTGGCTGCTGCCAGTGGCGCTTTGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAGCCGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTATGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAACGGCGCAGCCTTTTCCTGGTTCTCGTTTTTTGCTCACATGTTTCTTTTGGCGTTATCCCCTGATTCTGTGGATAACCGCATCTCCGCTTTTGAGTGAGCAGACACCGCTCGCCGCAGCCGAACGACCGAGTGTAGCGAGTCAGTGAGCGAGGAAGCGGAAGAGCGCCGGAACGTGCATTTTCTCCTTACGCATCTGTGCGGCATTTCACATCGGACATGGTGCGCTTTCCATACAATTCGTACTGATGCCGCATAGTTAAGCCAGTATACACTCCGCTATCGCTACGTGACTGGTTCAGGGCTTCGCCCCGAAACCCCCTGACGCGCCCTGAGGGGCTTGTCTGCTCCCGGCATCCGCTCACAGACAAGCTGTTACCGTCTCCGGGAGCTGTATGTGTCAGAGGTTTTCACCGTCATCCCCGAAGCGTGCGA"
population.fitness(seq)
population.select_bests()
population.do_gen()

def test_get_pop():
    li = population.get_pop()
    assert len(li) == 20
    assert type(li[0]) is Chromosome
test_get_pop()

def test_get_distance():
    res = population.get_distance()[0]
    assert type(res) is numpy.float64 or type(res) is float
test_get_distance()

def test_fitness():
    seq="AAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCCAGTAAACGAAAAAACCGCCTGGGGAGGCGGTTTAGTCGAAGGTTAAGTCAGTTGGGGACTGCTTAACCGGGTAACTGGCTTGGTGGAGCACAGATACCAAATACTGTCCTTCTAGTGTAGCCGCAGTTAGGCCACCACTTCAAGAACTCTTAATATCTCAATCCACCTTGTCCAGTTACCAGTGGCTGCTGCCAGTGGCGCTTTGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAGCCGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTATGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAACGGCGCAGCCTTTTCCTGGTTCTCGTTTTTTGCTCACATGTTTCTTTTGGCGTTATCCCCTGATTCTGTGGATAACCGCATCTCCGCTTTTGAGTGAGCAGACACCGCTCGCCGCAGCCGAACGACCGAGTGTAGCGAGTCAGTGAGCGAGGAAGCGGAAGAGCGCCGGAACGTGCATTTTCTCCTTACGCATCTGTGCGGCATTTCACATCGGACATGGTGCGCTTTCCATACAATTCGTACTGATGCCGCATAGTTAAGCCAGTATACACTCCGCTATCGCTACGTGACTGGTTCAGGGCTTCGCCCCGAAACCCCCTGACGCGCCCTGAGGGGCTTGTCTGCTCCCGGCATCCGCTCACAGACAAGCTGTTACCGTCTCCGGGAGCTGTATGTGTCAGAGGTTTTCACCGTCATCCCCGAAGCGTGCGA"
    population=Population(20)
    population.fitness(seq)
    res = population.get_distance()
    assert len(res) == 20
    assert type(res[0]) is numpy.float64
test_fitness()

def test_bests():
    population=Population(20)
    seq="AAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCCAGTAAACGAAAAAACCGCCTGGGGAGGCGGTTTAGTCGAAGGTTAAGTCAGTTGGGGACTGCTTAACCGGGTAACTGGCTTGGTGGAGCACAGATACCAAATACTGTCCTTCTAGTGTAGCCGCAGTTAGGCCACCACTTCAAGAACTCTTAATATCTCAATCCACCTTGTCCAGTTACCAGTGGCTGCTGCCAGTGGCGCTTTGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAGCCGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTATGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAACGGCGCAGCCTTTTCCTGGTTCTCGTTTTTTGCTCACATGTTTCTTTTGGCGTTATCCCCTGATTCTGTGGATAACCGCATCTCCGCTTTTGAGTGAGCAGACACCGCTCGCCGCAGCCGAACGACCGAGTGTAGCGAGTCAGTGAGCGAGGAAGCGGAAGAGCGCCGGAACGTGCATTTTCTCCTTACGCATCTGTGCGGCATTTCACATCGGACATGGTGCGCTTTCCATACAATTCGTACTGATGCCGCATAGTTAAGCCAGTATACACTCCGCTATCGCTACGTGACTGGTTCAGGGCTTCGCCCCGAAACCCCCTGACGCGCCCTGAGGGGCTTGTCTGCTCCCGGCATCCGCTCACAGACAAGCTGTTACCGTCTCCGGGAGCTGTATGTGTCAGAGGTTTTCACCGTCATCCCCGAAGCGTGCGA"
    population=Population(20)
    population.fitness(seq)
    population.select_bests()
    res = population.get_bests()
    assert len(res) == 20//10
    assert type(res[0].get_chr()) is dict
test_bests()

def test_do_gen():
    population=Population(20)
    seq="AAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCCAGTAAACGAAAAAACCGCCTGGGGAGGCGGTTTAGTCGAAGGTTAAGTCAGTTGGGGACTGCTTAACCGGGTAACTGGCTTGGTGGAGCACAGATACCAAATACTGTCCTTCTAGTGTAGCCGCAGTTAGGCCACCACTTCAAGAACTCTTAATATCTCAATCCACCTTGTCCAGTTACCAGTGGCTGCTGCCAGTGGCGCTTTGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAGCCGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTATGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAACGGCGCAGCCTTTTCCTGGTTCTCGTTTTTTGCTCACATGTTTCTTTTGGCGTTATCCCCTGATTCTGTGGATAACCGCATCTCCGCTTTTGAGTGAGCAGACACCGCTCGCCGCAGCCGAACGACCGAGTGTAGCGAGTCAGTGAGCGAGGAAGCGGAAGAGCGCCGGAACGTGCATTTTCTCCTTACGCATCTGTGCGGCATTTCACATCGGACATGGTGCGCTTTCCATACAATTCGTACTGATGCCGCATAGTTAAGCCAGTATACACTCCGCTATCGCTACGTGACTGGTTCAGGGCTTCGCCCCGAAACCCCCTGACGCGCCCTGAGGGGCTTGTCTGCTCCCGGCATCCGCTCACAGACAAGCTGTTACCGTCTCCGGGAGCTGTATGTGTCAGAGGTTTTCACCGTCATCCCCGAAGCGTGCGA"
    population=Population(20)
    population.fitness(seq)
    population.select_bests()
    population.do_gen()
    res = population.do_gen()
    assert type(res) is Population
test_do_gen()

def test_crossover():
    population=Population(20)
    seq="AAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCCAGTAAACGAAAAAACCGCCTGGGGAGGCGGTTTAGTCGAAGGTTAAGTCAGTTGGGGACTGCTTAACCGGGTAACTGGCTTGGTGGAGCACAGATACCAAATACTGTCCTTCTAGTGTAGCCGCAGTTAGGCCACCACTTCAAGAACTCTTAATATCTCAATCCACCTTGTCCAGTTACCAGTGGCTGCTGCCAGTGGCGCTTTGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAGCCGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTATGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAACGGCGCAGCCTTTTCCTGGTTCTCGTTTTTTGCTCACATGTTTCTTTTGGCGTTATCCCCTGATTCTGTGGATAACCGCATCTCCGCTTTTGAGTGAGCAGACACCGCTCGCCGCAGCCGAACGACCGAGTGTAGCGAGTCAGTGAGCGAGGAAGCGGAAGAGCGCCGGAACGTGCATTTTCTCCTTACGCATCTGTGCGGCATTTCACATCGGACATGGTGCGCTTTCCATACAATTCGTACTGATGCCGCATAGTTAAGCCAGTATACACTCCGCTATCGCTACGTGACTGGTTCAGGGCTTCGCCCCGAAACCCCCTGACGCGCCCTGAGGGGCTTGTCTGCTCCCGGCATCCGCTCACAGACAAGCTGTTACCGTCTCCGGGAGCTGTATGTGTCAGAGGTTTTCACCGTCATCCCCGAAGCGTGCGA"
    population=Population(20)
    population.fitness(seq)
    population.select_bests()
    population.do_gen()
    aux=Chromosome()
    res = aux.crossover(population.get_indiv(0),population.get_indiv(1))
    assert type(res) is tuple
test_crossover()

def test_mutation():
    population=Population(20)
    print(population.get_pop())
    seq="AAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCCAGTAAACGAAAAAACCGCCTGGGGAGGCGGTTTAGTCGAAGGTTAAGTCAGTTGGGGACTGCTTAACCGGGTAACTGGCTTGGTGGAGCACAGATACCAAATACTGTCCTTCTAGTGTAGCCGCAGTTAGGCCACCACTTCAAGAACTCTTAATATCTCAATCCACCTTGTCCAGTTACCAGTGGCTGCTGCCAGTGGCGCTTTGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAGCCGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTATGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAACGGCGCAGCCTTTTCCTGGTTCTCGTTTTTTGCTCACATGTTTCTTTTGGCGTTATCCCCTGATTCTGTGGATAACCGCATCTCCGCTTTTGAGTGAGCAGACACCGCTCGCCGCAGCCGAACGACCGAGTGTAGCGAGTCAGTGAGCGAGGAAGCGGAAGAGCGCCGGAACGTGCATTTTCTCCTTACGCATCTGTGCGGCATTTCACATCGGACATGGTGCGCTTTCCATACAATTCGTACTGATGCCGCATAGTTAAGCCAGTATACACTCCGCTATCGCTACGTGACTGGTTCAGGGCTTCGCCCCGAAACCCCCTGACGCGCCCTGAGGGGCTTGTCTGCTCCCGGCATCCGCTCACAGACAAGCTGTTACCGTCTCCGGGAGCTGTATGTGTCAGAGGTTTTCACCGTCATCCCCGAAGCGTGCGA"
    population=Population(20)
    population.fitness(seq)
    population.select_bests()
    population.do_gen()
    aux=Chromosome()
    aux.apply_mutation()
    res = population.get_pop()
    assert type(res[0]) is Chromosome and type(aux) is Chromosome
test_mutation()

