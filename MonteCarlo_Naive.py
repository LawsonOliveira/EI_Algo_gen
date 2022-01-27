from MonteCarlo import *
import numpy as np
from math import sqrt
from RotTable import RotTable
from numpy import angle
# from EI_Algo_gen.RotTable import RotTable
import time
root=node(node.__ORIGINAL_ROT_TABLE)

filename='plasmid_8k.fasta'
lineList = [line.rstrip('\n') for line in open(filename)]
seq = ''.join(lineList[1:])

def createrandomtable():

