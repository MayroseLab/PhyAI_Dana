import re, os, shutil
import pandas as pd
from ete3 import Tree
from subprocess import Popen, PIPE, STDOUT

RAXML_NG_SCRIPT = "raxml-ng"    # after you install raxml-ng on your machine
MSA_PHYLIP_FILENAME = "masked_species_real_msa.phy"



def prune_branch(t_orig, prune_name):
	'''
	returns (a copy of) both ETE subtrees after pruning
	'''
	t_cp_p = t_orig.copy()  # the original tree is needed for each iteration
	assert t_cp_p & prune_name    # todo Oz: add indicative error
	prune_node_cp = t_cp_p & prune_name  # locate the node in the copied subtree
	assert prune_node_cp.up

	nname = prune_node_cp.up.name
	prune_loc = prune_node_cp
	prune_loc.detach()  # pruning: prune_node_cp is now the subtree we detached. t_cp_p is the one that was left behind
	t_cp_p.search_nodes(name=nname)[0].delete(preserve_branch_length=True)  # delete the specific node (without its childs) since after pruning this branch should not be divided
	return nname, prune_node_cp, t_cp_p


def regraft_branch(t_cp_p, prune_node_cp, rgft_name, nname, preserve=False):
	'''
	recieves: 2 ETE subtrees and 2 node names
	returns: an ETEtree with the 2 concatenated ETE subtrees
	'''
	t_temp = Tree()  # for concatenation of both subtrees ahead, to avoid polytomy
	t_temp.add_child(prune_node_cp)
	t_curr = t_cp_p.copy()
	rgft_node_cp = t_curr & rgft_name  # locate the node in the copied subtree
	new_branch_length = rgft_node_cp.dist / 2

	rgft_loc = rgft_node_cp.up
	rgft_node_cp.detach()
	t_temp.add_child(rgft_node_cp, dist=new_branch_length)
	t_temp.name = nname
	rgft_loc.add_child(t_temp, dist=new_branch_length)  # regrafting
	if nname == "ROOT_LIKE":  # (4)
		t_temp.delete()
		# preserve the name of the root node, as this is a REAL node in this case
		preserve = True
	return t_curr, preserve



def SPR_generator(t):
	'''
	(1) avoid moves to a branch in the pruned subtree
	(2) avoid moves to the branch of a sibiling (retaining the same topology)
	(3) avoid moves to the branch leading to the parent of the (retaining the same topology)
	(4) handle the automatic "ROOTLIKE" node of ETE trees -
	    if the PRUNE location is around the rootlike --> delete the "ROOT_LIKE" subnode when pasting in dest,
	    and preserve the name of the REAL rootnode when converting back to newick
	'''
	
	prune_name = "N8"
	rgft_name = "N4"
	all_nodes_names = [n.name for n in t.get_descendants()]
	
	# perform pruning
	nname, subtree1, subtree2 = prune_branch(t, prune_name)  # subtree1 is the pruned subtree. subtree2 is the remaining subtree
	if nname == (t&rgft_name).up.name:  # captures (2)
		pass
		#continue
	remaining_nodes_names = [n.name for n in subtree2.get_descendants()]
	if rgft_name in remaining_nodes_names:  # captures (1) and (3)
		neighbor_tree, preserve = regraft_branch(subtree2, subtree1, rgft_name, nname)
		if preserve:  # namely, (4) is True
			neighbor_tree_str = neighbor_tree.write(format=1, format_root_node=True)
		else:
			neighbor_tree_str = neighbor_tree.write(format=1)
		print(neighbor_tree_str)
		print(Tree(newick=neighbor_tree_str, format=1).get_ascii(show_internal=True))
	else:
		print("Forbidden move: {} --> {} for this tree !".format(prune_name, rgft_name))
		#continue
		
	return





if __name__ == '__main__':
	from ete3 import Tree, TreeStyle, TextFace
	
	t = Tree("((a,b),c);")
	# Basic tree style
	ts = TreeStyle()
	# ts.show_leaf_name = True
	
	# Add two text faces to different columns
	t.add_face(TextFace("hola "), column=0, position="branch-right")
	t.add_face(TextFace("mundo!"), column=1, position="branch-right")
	t.show(tree_style=ts)
	exit()
	
	tree_str_with_internal_names = "(((((Sp002:2.9e-07,Sp006:0.00173012)N8:0.00248844,Sp005:0.00978979)N6:0.0136241,Sp004:0.0388109)N4:0.0316201,Sp000:0.045075)N1:0.0837587,Sp003:0.00308638,Sp001:0.0180697);"
	t = Tree(newick=tree_str_with_internal_names, format=1)
	t.get_tree_root().name = "ROOT_LIKE"
	print(t.get_ascii(show_internal=True))
	SPR_generator(t)