from Algo_gen.Genetic_Algo import *
import argparse
from os import popen
from re import A
from Algo_MCTS.MCTS import *

parser = argparse.ArgumentParser()
parser.add_argument("--filename", help="input filename of DNA sequence")
parser.parse_args()
args = parser.parse_args()

# Command pour running with a file : python main.py --filename=plasmid_8k.fasta

def main(p=2):
    # Running of genetics algorithm
    seq="AAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCCAGTAAACGAAAAAACCGCCTGGGGAGGCGGTTTAGTCGAAGGTTAAGTCAGTTGGGGACTGCTTAACCGGGTAACTGGCTTGGTGGAGCACAGATACCAAATACTGTCCTTCTAGTGTAGCCGCAGTTAGGCCACCACTTCAAGAACTCTTAATATCTCAATCCACCTTGTCCAGTTACCAGTGGCTGCTGCCAGTGGCGCTTTGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAGCCGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTATGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAACGGCGCAGCCTTTTCCTGGTTCTCGTTTTTTGCTCACATGTTTCTTTTGGCGTTATCCCCTGATTCTGTGGATAACCGCATCTCCGCTTTTGAGTGAGCAGACACCGCTCGCCGCAGCCGAACGACCGAGTGTAGCGAGTCAGTGAGCGAGGAAGCGGAAGAGCGCCGGAACGTGCATTTTCTCCTTACGCATCTGTGCGGCATTTCACATCGGACATGGTGCGCTTTCCATACAATTCGTACTGATGCCGCATAGTTAAGCCAGTATACACTCCGCTATCGCTACGTGACTGGTTCAGGGCTTCGCCCCGAAACCCCCTGACGCGCCCTGAGGGGCTTGTCTGCTCCCGGCATCCGCTCACAGACAAGCTGTTACCGTCTCCGGGAGCTGTATGTGTCAGAGGTTTTCACCGTCATCCCCGAAGCGTGCGA"
    if p==1:
        pop_size=40  # It has to be bigger than 20 and a multiple of 4, for example, 40
        algo=Genetic_Algo(pop_size,100,seq,0.01,args)
        algo.darwin()
    # Running the Monte Carlo algorithm
    else: 
        a=MCTS(500,seq)
        a.draw()
        
if __name__ == "__main__" :
    main(1)
