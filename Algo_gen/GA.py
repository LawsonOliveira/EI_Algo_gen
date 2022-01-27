from os import popen
from re import A
from Population import *

class GA:
    def __init__(self,n,max_iter,seq,dist_max,args=None):
        self.__pop_size = n
        self.__max_iter=max_iter
        self.__args=args
        self.__max_dist=dist_max
        if self.__args.filename:
            # Read file
            lineList = [line.rstrip('\n') for line in open(self.__args.filename)]
            # Formatting
            self.__seq = ''.join(lineList[1:])
        else:
            self.__seq=seq


    def darwin(self):
        new_distance=math.inf
        p=0
        population=Population(self.__pop_size)
        while new_distance>self.__max_dist and p<self.__max_iter:
            p+=1
            population.select_bests()
            population.do_gen()
            population.mutation()
            population.fitness(self.__seq)
            best,new_distance=population.pickbest()
            print("Iteration",p)
            print("Distance :",new_distance)


        traj=Traj3D()
        if self.__args.filename:
            archive = open(self.__args.filename+"_solution_.fasta",'w')
            archive.write(str(best.get_chr()))
            archive.close()
            traj.compute(self.__seq,best)
            traj.draw(self.__args.filename+".png")

        else:
            archive = open("exemple"+"_solution_.fasta",'w')
            archive.write(str(best.get_chr()))
            archive.close()
            traj.compute(self.__seq,best)
            traj.draw("sample.png")

