import sys  
from pathlib import Path  
file = Path(__file__).resolve()  
package_root_directory = file.parents[1]  
sys.path.append(str(package_root_directory))

from Algo_MCTS.Monte_Carlo_fonctions import *
from Algo_MCTS.Traj3D import *
from Algo_MCTS.Node import *

# Main classe for the Monte Carlo Tree Search
class MCTS:
    
    def __init__(self, nbit,seq,critere=10**-3):
        self.best = MCTSfun(nbit, seq,critere)
        self.seq=seq
        self.traj=Traj3D()

    def draw(self):
        a = Node()
        a.newTable(self.best)
        self.traj.compute(self.seq,a)
        self.traj.draw("MTRS.png")

    def getsol(self):
        return self.best  # ICI LOUCHE
    
    def getTraj(self):
        return self.traj