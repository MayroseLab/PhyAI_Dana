import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")

from defs import *
from collections import Counter
from utils.msa_functions import *
from itertools import combinations

#from os.path import exists
#from ete3 import Tree
#from collections import OrderedDict



def cp_internal_names(treepath_no_internal, treepath_with_internal):
	with open(treepath_with_internal) as fp:
		with_names = fp.read()
	with open(treepath_no_internal) as fp:
		nonames = fp.read()

	with_names_bls = re.findall(":(\d*\.?\d+)", with_names)
	nonames_bls = re.findall(":(\d*\.?\d+)", nonames)
	while len(set(with_names_bls)) != len(with_names_bls):
		u = [k for (k, v) in Counter(with_names_bls).items() if v > 1][0]
		ix = with_names_bls.index(u)
		with_names_bls[ix] = u + "1"
		with_names = with_names.replace(u, u + "1", 1)


	dict = {with_names_bls[i]: nonames_bls[i] for i in range(len(with_names_bls))}

	regex = re.compile("(%s)" % "|".join(map(re.escape, dict.keys())))
	try:
		new_str = regex.sub(lambda mo: dict[mo.string[mo.start():mo.end()]], with_names)
	except:
		print("in func cp_internal_names:")
		print(treepath_no_internal)
		print(dict)

	with open(treepath_no_internal, 'r') as fp:
		tree_str = fp.read()
		if re.search("\)N", tree_str):
			return
	with open(treepath_no_internal, 'w') as fp:
		fp.write(new_str)


def get_newick_tree(tree):
	"""
	:param tree: newick tree string or txt file containing one tree
	:return:	tree: a string of the tree in ete3.Tree format
	"""
	if os.path.exists(tree):
		with open(tree, 'r') as tree_fpr:
			tree = tree_fpr.read().strip()
	return tree


def get_branch_lengths(tree):
	"""
	:param tree: Tree node or tree file or newick tree string;
	:return: list of branch lengths
	"""
	try:
		if type(tree) == str:
			tree = Tree(get_newick_tree(tree), format=1)
		tree_root = tree.get_tree_root()
	except:
		print(tree)
	if len(tree) == 1 and not "(" in tree:  # in one-branch trees, sometimes the newick string is without "(" and ")" so the .iter_decendants returns None
		return [tree.dist]
	branches = []
	for node in tree_root.iter_descendants(): # the root dist is 1.0, we don't want it
		branches.append(node.dist)

	return branches


def get_total_branch_lengths(tree):
	"""
	:param tree: Tree node or tree file or newick tree string;
	:return: total branch lengths
	"""
	branches = get_branch_lengths(tree)
	return sum(branches)
	

def estimate_lengths(t, rgft_node):
	res_tbl = get_total_branch_lengths(t)   # tbl after naive bl estimations (like in RAxML) in both the pruning and the rgft locations. before uptimization!
	res_bl = rgft_node.dist   				# which is half of the orig rgft branch length
	return res_tbl, res_bl
	

def dist_between_nodes(t, node1):
	nleaves_between, tbl_between = {},{}
	for node2 in t.get_descendants("levelorder")[::-1]:
		if not node2.name:
			node2.name = "Nnew"
		nname2 = node2.name
		if node1.name == nname2:
			continue

		nleaves_between[nname2] = node1.get_distance(node2, topology_only=True)+1   # +1 to convert between nodes count to edges
		tbl_between[nname2] = node1.get_distance(node2, topology_only=False)

	return nleaves_between, tbl_between


def init_recursive_features(t):
	assert isinstance(t, Tree)
	for node in t.traverse("postorder"):
		if node.is_leaf():
			node.add_feature("cumBL", 0)
			node.add_feature("maxBL", 0)
			node.add_feature("ntaxa", 1)
		else:
			#since it's postorder, we must have already went over its children leaves
			left, right = node.children
			node.add_feature("cumBL", left.cumBL + right.cumBL + left.dist + right.dist)
			node.add_feature("maxBL", max(left.maxBL, right.maxBL, left.dist, right.dist))
			node.add_feature("ntaxa", left.ntaxa + right.ntaxa)


def update_node_features(subtree):
	"""
	:param subtree: a node that needs update. might be None or a leaf
	:return: None
	"""
	# if subtree and not subtree.is_leaf():
	left, right = subtree.children
	subtree.cumBL = left.cumBL + right.cumBL + left.dist + right.dist
	subtree.maxBL = max(left.maxBL, right.maxBL, left.dist, right.dist)
	subtree.ntaxa = left.ntaxa + right.ntaxa

"""
def dist_between_nodes_mat(t):
	'''
	:param tree obj:
	:return: distances between nodes matrix
	'''
	import pandas as pd
	df = pd.DataFrame()

	


	for i in range(ntaxa-1):
		df.loc[msa[i].id, msa[i].id] = 0
		for j in range(i+1, ntaxa):
			pdist = hamming_distance(msa[i], msa[j])
			df.loc[msa[i].id, msa[j].id] = pdist/nchars
			df.loc[msa[j].id, msa[i].id] = pdist/nchars
	df.loc[msa[-1].id, msa[-1].id] = 0

	return nleaves_between, tbl_between
"""

def dist_between_nodes(t, node1):
	### todo: change! make an NxN matrix for all nodes of the complete tree (outide) and compute the distances between each pair
	nleaves_between, tbl_between = {},{}
	for node2 in t.get_descendants("levelorder")[::-1]:
		if not node2.name:
			node2.name = "Nnew"
		nname2 = node2.name
		if node1.name == nname2:
			continue

		nleaves_between[nname2] = node1.get_distance(node2, topology_only=True)+1   # +1 to convert between nodes count to edges
		tbl_between[nname2] = node1.get_distance(node2, topology_only=False)

	return nleaves_between, tbl_between


def prune_branch(t_orig, prune_name):
	'''
	get (a copy of) both subtrees after pruning
	'''
	t_cp_p = t_orig.copy()  				# the original tree is needed for each iteration
	prune_node_cp = t_cp_p & prune_name     # locate the node in the copied subtree
	assert prune_node_cp.up

	nname = prune_node_cp.up.name
	prune_loc = prune_node_cp
	prune_loc.detach()  # pruning: prune_node_cp is now the subtree we detached. t_cp_p is the one that was left behind
	t_cp_p.search_nodes(name=nname)[0].delete(preserve_branch_length=True)  # delete the specific node (without its childs) since after pruning this branch should not be divided
	if not nname:
		nname = "Nnew_p"

	return nname, prune_node_cp, t_cp_p



def calc_leaves_features(tree_str, move_type, rgft_node_name=None):
	
	if not move_type == "res":
		t = Tree(newick=tree_str, format=1)
	
		name2bl, name2pdist_pruned, name2pdist_remaining, name2tbl_pruned, name2tbl_remaining, name2longest_pruned, name2longest_remaining, name2ntaxa,name2ntaxa_pruned, name2ntaxa_remaining, name2pars_pruned, name2parse_remaining,names2topology_dist, names2bl_dist = {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}
		for node in t.get_descendants("levelorder")[::-1]:
			nname = node.name
			namex, subtree1, subtree2 = prune_branch(t, nname)
	
			name2bl[nname] = node.dist
	
			name2tbl_pruned[nname] = get_total_branch_lengths(subtree1)
			name2tbl_remaining[nname] = get_total_branch_lengths(subtree2)
			name2longest_pruned[nname] = max(get_branch_lengths(subtree1))
			name2longest_remaining[nname] = max(get_branch_lengths(subtree2))
			name2ntaxa_pruned[nname] = len(subtree1.get_tree_root())
			name2ntaxa_remaining[nname] = len(subtree2.get_tree_root())
			
			res_dists = ["unrelevant", "unrelevant"] if move_type == "rgft" else dist_between_nodes(t, node)
			names2topology_dist[nname] = res_dists[0]	# res_dists is a dict
			names2bl_dist[nname] = res_dists[1]			# res_dists is a dict
			
		d = OrderedDict([("bl", name2bl), ("longest", max(get_branch_lengths(t))),
		                 ("ntaxa_p", name2ntaxa_pruned), ("ntaxa_r", name2ntaxa_remaining),
		                 ("tbl_p", name2tbl_pruned), ("tbl_r", name2tbl_remaining),
						 ("top_dist", names2topology_dist), ("bl_dist", names2bl_dist),
		                 ("longest_p", name2longest_pruned), ("longest_r", name2longest_remaining),])


	else:
		t = Tree(tree_str, format=1)  # the res tree NOT optimized
		rgft_node = t & rgft_node_name
		res_tbl, res_bl = estimate_lengths(t, rgft_node)
		d = OrderedDict([("res_tbl", res_tbl), ("res_bl", res_bl)])
	
	##print (time_pdist, time_bl, time_tbl, time_longest, time_ntaxa, time_pars, time_dists, time_res_bls)
	#data_res = pd.DataFrame.from_dict(d)
	#pd.set_option('display.max_columns', 50)
	#print(data_res)
	##data_res.to_csv(DATA_PATH + "blah.csv")
	##print(d)
	return d


def calc_leaves_features_Shiran(tree_str, move_type, rgft_node_name=None):
	if not move_type == "res":
		t = Tree(newick=tree_str, format=1)
		name2bl, name2pdist_pruned, name2pdist_remaining, name2tbl_pruned, name2tbl_remaining, name2longest_pruned, name2longest_remaining, name2ntaxa, name2ntaxa_pruned, name2ntaxa_remaining, name2pars_pruned, name2parse_remaining, names2topology_dist, names2bl_dist = {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}

		# first get tree leaves+nodes - and preserve the order. Do not unroot at any point.
		# set one leaf for a primary rooting
		nodes_order = []
		first_node = t.get_leaves()[0]
		#print(first_node.name)
		t.set_outgroup(first_node)


		# We want to traverse the tree by pooling the parent of the remaining tree and moving it to the outgroup -
		# so we do it in a preorder manner, such that
		for node in t.traverse("preorder"):
			nodes_order.append(node)

		# initialize tree features in here, and afterwards only rotate tree and update them
		nodes_order.pop(0)  # this is t, first_node's nonsense parent
		init_recursive_features(t)  # compute features for primary tree
		for node in nodes_order:
			#if node is first_node:    # EDITED DANA
			#	continue              # EDITED DANA
			#else:
			# root the tree with current node's parent
			# when rotating a tree in a preorder manner - the parent's parent (grandparent) becomes the parent's son.
			# so we need to update both of them (if exist)
			# print(t.get_ascii(show_internal=True))
			nodes_to_update = [node.up]
			while nodes_to_update[-1]:
				# print(nodes_to_update[-1].name)
				nodes_to_update.append(nodes_to_update[-1].up)
			nodes_to_update.pop(-1)  # None
			nodes_to_update.pop(-1)  # the nonesense root

			t.set_outgroup(node)
			for up_node in nodes_to_update[::-1]:
				update_node_features(up_node)

			nname = node.name
			print(nname)
			if not "Sp" in nname:                      # EDITED DANA
				subtree1, subtree2 = node.children     # EDITED DANA
			else:
				subtree1, subtree2= node, node
			print(node.get_ascii(show_internal=True))
			print(subtree1.get_ascii(show_internal=True))
			print(subtree2.get_ascii(show_internal=True))
			print("######3")
			name2bl[nname] = subtree2.dist * 2 if "Sp" in nname else subtree2.dist*2  # EDITED DANA
			name2tbl_pruned[nname] = subtree1.cumBL  + subtree1.dist * 2 # EDITED DANA
			name2tbl_remaining[nname] = subtree2.cumBL
			name2longest_pruned[nname] = max(subtree1.maxBL, subtree1.dist * 2)
			name2longest_remaining[nname] = subtree2.maxBL
			name2ntaxa_pruned[nname] = subtree1.ntaxa
			name2ntaxa_remaining[nname] = subtree2.ntaxa

		# todo: compute once for the complete tree in the external function
		# res_dists = ["unrelevant", "unrelevant"] if move_type == "rgft" else dist_between_nodes(t, node)
		# names2topology_dist[nname] = res_dists[0]  # res_dists is a dict
		# names2bl_dist[nname] = res_dists[1]  # res_dists is a dict

		d = OrderedDict([("bl", name2bl), ("longest", max(subtree1.maxBL, subtree2.maxBL, subtree2.dist * 2)),
						 ("ntaxa_p", name2ntaxa_pruned), ("ntaxa_r", name2ntaxa_remaining),
						 ("tbl_p", name2tbl_pruned), ("tbl_r", name2tbl_remaining),
						 # ("top_dist", names2topology_dist), ("bl_dist", names2bl_dist),
						 ("longest_p", name2longest_pruned), ("longest_r", name2longest_remaining)])


	else:
		t = Tree(tree_str, format=1)
		rgft_node = t & rgft_node_name
		res_tbl, res_bl = estimate_lengths(t, rgft_node)
		d = OrderedDict([("res_tbl", res_tbl), ("res_bl", res_bl)])

	data_res = pd.DataFrame.from_dict(d)
	pd.set_option('display.max_columns', 50)
	print(data_res)
	return d



if __name__ == '__main__':
	#calc_leaves_features(DIRPATH + PHYML_TREE_FILENAME.format('bionj'), DIRPATH + MSA_PHYLIP_FILENAME)
	pass