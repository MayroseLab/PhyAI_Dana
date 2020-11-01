import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")


from defs import *

from subprocess import Popen, PIPE, STDOUT
RAXML_NG_SCRIPT = "raxml-ng"












if __name__ == '__main__':
	pass
	t = Tree("/groups/itay_mayrose/danaazouri/PhyAI/starting_trees_ml_minus1/data/training_datasets/Sp0019_Sp0023/bad_tree.tre", format=1)
	t.set_outgroup(t&"Sp0025")
	(t&"Sp0025").unroot()

	print(t.write(format=1))
