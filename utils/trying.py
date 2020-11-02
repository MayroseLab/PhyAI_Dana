import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")


from defs import *

from subprocess import Popen, PIPE, STDOUT
RAXML_NG_SCRIPT = "raxml-ng"












if __name__ == '__main__':
	pass
	t = Tree("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/data/training_datasets/exampleSphaero/masked_species_real_msa.phy_phyml_tree_bionj.txt", format=1)
	print(t.get_ascii(show_internal=True))