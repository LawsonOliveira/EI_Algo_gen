# ST2 (Game Theory) - IE Genetic Algorithms

# Description 

Created two programs to make um plasmid circular using object oriented programming.

# Team members: 
- Lawson OLIVEIRA LIMA
- Ahmed EL BAJDALI
- Pierrick BOURNEZ 
- osvaldo CARTAGENA

# Project: Coding a sudoku game

Project goal: The objective is to use two different approaches to render a circular DNA sequence, more specifically, the Monte Carlo Tree Search algorithm and genetic algorithms.

###########################################################################################################

# Here are the modules that need to be imported to use our project

- mathutils
- numpy
- matplotlib

# Using our project

To use our project, you need to run the main file and choose the algorithm you want to use, MCTS or Genetic Algorithm

###########################################################################################################
# Organization of files and modules

- Algo_gen : this file contains the codes of the Genetic Algorithm
    - Benchmark_GA : this file contains our benchmark class for the genetic ago
    - Chromosome : this file contains our Chromosoome, it represents a set of 16 genes
    - Population : this file contains our Population, it represents a set of 100 chromosomes
    - Genetics Algo : this file contains our main class used by the other classes, it makes the link between the main functions
    - Traj3D : this file contains our function that applies each individual of the population (our chromosome) in the dna sequence
    - Test : this file contains all the tests we made for Genetic Algorithm

- Solutions_Genetics_Algo : this folder contains the solutions obtained using a genetic algorithm

- Solutions_Monte_Carlo : this folder contains the solutions obtained using the MCTS

- Algo_MCTS : this folder contains the MCTS codes
    - Benchmark_MCTS : this file contains our benchmark class for the Monte Carlo Tree Search
    - Monte_Carlo_functions : this file contains the functions used for the MCTSomosoome, it represents a set of 16 genes
    - MCTS : this file contains our main class used by the other classes, it makes the link between the main functions
    - Node : this file contains our class that defines the nodes of the tree

- Algo_MCTS : this file contains the function test


###########################################################################################################

# Our strong points
    - Two very powerful programs that obtain optimal solutions
    - An efficient way to do the crossover 
    - Objective functions that increase the quality of the solution


###########################################################################################################

# Suggestions to improve our project in the future
    - Slightly improve the fit function 
    - The constant c for the fit function should be more deeply analywe 
    - As the fit function is differentiable, some other algorithm would be more efficient 

# The results of the genetic algorithm were obtained for a population of 40 individuals
