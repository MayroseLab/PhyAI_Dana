import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")
import warnings

from defs import *
from Bio import Phylo
import networkx, pylab
TREE_PATH = "/groups/itay_mayrose/danaazouri/PhyAI/DBset2/data_test/tree.txt"

tree = Phylo.read(TREE_PATH, "newick")
net = Phylo.to_networkx(tree)
pos = networkx.spring_layout(net)
print(pos)


networkx.draw(net)
pylab.show()