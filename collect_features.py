#########################################################################
##                 Copyright (C). All Rights Reserved.                   ##
##      "Harnessing machine-learning to boost heuristic strategies       ##
##                                      for phylogenetic-tree search"    ##
##                                                                       ##
## by Dana Azouri, Shiran Abadi, Yishay Mansour, Itay Mayrose, Tal Pupko ##
##                                                                       ##
##                                                                       ##
##   For information, contact danaazouri@mail.tau.ac.il                  ##
##                                                                       ##
## For academic, non-commercial use.                                     ##
## If you use the code, please cite the paper                            ##
##                                                                       ##
#########################################################################

from defs_PhyAI import *

################################################################################################
########################## begining of 'features extraction' section ###########################
################################################################################################

def get_msa_from_file(msa_file_path):
	#open file if exists
	if not os.path.exists(msa_file_path):
		return None
	try:
		msa = AlignIO.read(msa_file_path, PHYLIP_FORMAT)
	except:
		return None
	return msa


def get_msa_properties(msa):
	"""
	:param msa: bio.AlignIO format or path to msa file
	:return:
	"""
	if isinstance(msa, str):
		msa = get_msa_from_file(msa)
	ntaxa = len(msa)
	nchars = msa.get_alignment_length()

	return ntaxa, nchars

	
def get_seqs_dict(msa_file):
	alignment = get_msa_from_file(msa_file)
	seqs_dict = {seq.id: seq.seq for seq in alignment}
	return seqs_dict


def trunc_msa(subt, seqs_dict, trunc_msa_path=None):
	if type(subt) == str:
		subt = Tree(newick=subt, format=1)
	records = []
	for leaf in subt.iter_leaves():  # go over all leaves (all species) in this subtree
		leaf_name = leaf.name
		records.append(SeqRecord(seqs_dict[leaf_name], id=leaf_name))
	trunc_msa = MultipleSeqAlignment(records)
	if trunc_msa_path:
		AlignIO.write(trunc_msa, trunc_msa_path, PHYLIP_FORMAT)
		rewrite_in_phylip(trunc_msa_path)

	return trunc_msa


def rewrite_in_phylip(msa_file):
	data = []
	with open(msa_file, 'r') as fp:
		for line in fp.readlines():
			nline = line.rstrip("\r\n")
			re_name_only = re.search("^(\S+)\s+\S+",nline)
			if re_name_only:
				name_only = re_name_only.group(1)
				end_name_ix = len(name_only) +1
				with_spaces = nline[:end_name_ix] + "      " + nline[end_name_ix:]
				nline = with_spaces
			data.append(nline)

	with open(msa_file, 'w') as nf:
		nf.write("\n".join(data))


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
	if len(tree) == 1 and not "(" in tree:    # in one-branch trees, sometimes the newick string is without "(" and ")" so the .iter_decendants returns None
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
	res_tbl = get_total_branch_lengths(t)  # tbl after naive bl estimations (like in RAxML) in both the pruning and the rgft locations. before uptimization!
	res_bl = rgft_node.dist				   # which is half of the orig rgft branch length
	
	return res_tbl, res_bl


def prune_branch(t_orig, prune_name): # the same function from 'SPR_and_lls' script
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



def calc_leaves_features(tree_str, move_type, rgft_node_name=None):
	if not move_type == "res":
		t = Tree(newick=tree_str, format=1)
		name2bl, name2pdist_pruned, name2pdist_remaining, name2tbl_pruned, name2tbl_remaining, name2longest_pruned, name2longest_remaining, name2ntaxa, name2ntaxa_pruned, name2ntaxa_remaining, name2pars_pruned, name2parse_remaining, names2topology_dist, names2bl_dist = {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}
		
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
			names2topology_dist[nname] = res_dists[0]  # res_dists is a dict
			names2bl_dist[nname] = res_dists[1]  # res_dists is a dict
		
		d = OrderedDict([("bl", name2bl), ("longest", max(get_branch_lengths(t))),
		                 ("ntaxa_p", name2ntaxa_pruned), ("ntaxa_r", name2ntaxa_remaining),
		                 ("tbl_p", name2tbl_pruned), ("tbl_r", name2tbl_remaining),
		                 ("longest_p", name2longest_pruned), ("longest_r", name2longest_remaining),
		                 ("top_dist", names2topology_dist), ("bl_dist", names2bl_dist)])
	
	else:
		t = Tree(newick=tree_str, format=1)  # res tree before optimization. replaced from des_gas
		rgft_node = t & rgft_node_name
		res_tbl, res_bl = estimate_lengths(t, rgft_node)
		d = OrderedDict([("res_tbl", res_tbl), ("res_bl", res_bl)])
	
	return d

################################################################################################
############################# end of 'features extraction' section #############################
################################################################################################
################################################################################################
##################### begining of 'organizing features in csv' section #########################
################################################################################################

def index_additional_rgft_features(df_rgft, ind, prune_name, rgft_name, features_restree_dict, features_dict_prune):
	df_rgft.ix[ind, FEATURES["top_dist"]] = features_dict_prune['top_dist'][prune_name][rgft_name]
	df_rgft.ix[ind, FEATURES["bl_dist"]] = features_dict_prune['bl_dist'][prune_name][rgft_name]
	df_rgft.ix[ind, FEATURES["res_bl"]] = features_restree_dict['res_bl']
	df_rgft.ix[ind, FEATURES["res_tbl"]] = features_restree_dict['res_tbl']

	return df_rgft


def index_shared_features(dff, ind, edge, move_type, features_dicts_dict):
	d_ll = format(dff.loc[ind, "ll"] - dff.loc[ind, "orig_ds_ll"], '.20f')
	dff.ix[ind, LABEL.format(move_type)] = d_ll  # LABEL

	#*tbl of orig tree will be calculated via 'RandomForest_learning' script		#f1  (FeaturePruning)
	dff.ix[ind, FEATURES["bl"]] = features_dicts_dict["bl"][edge] 					#f2
	dff.ix[ind, FEATURES["longest"]] = features_dicts_dict["longest"]				#f3


	#index subtrees features
	for subtype in ["p", "r"]:
		dff.ix[ind, FEATURES["ntaxa_{}".format(subtype)]] = features_dicts_dict["ntaxa_{}".format(subtype)][edge]  		#f4,5
		dff.ix[ind, FEATURES["tbl_{}".format(subtype)]] = features_dicts_dict["tbl_{}".format(subtype)][edge]  			#f6,7
		dff.ix[ind, FEATURES["longest_{}".format(subtype)]] = features_dicts_dict["longest_{}".format(subtype)][edge]	#f8,9

	return dff


def collect_features(ds_path, step_number, outpath_prune, outpath_rgft, tree_type='bionj'):
	dfr = pd.read_csv(TREES_PER_DS.format(ds_path, step_number), index_col=0)
	df_prune = pd.read_csv(outpath_prune, index_col=0)
	df_rgft = pd.read_csv(outpath_rgft, index_col=0)

	features_prune_dicts_dict = calc_leaves_features(ds_path + PHYML_TREE_FILENAME.format(tree_type), "prune")

	for i, row in dfr.iterrows():
		ind = row.name
		print(ind)
		tree = row["newick"]
		if row["rgft_name"] == SUBTREE2:	# namely the remaining subtree
			features_rgft_dicts_dict = calc_leaves_features(tree, "rgft")  # msa will be truncated within the function
		if not "subtree" in row["rgft_name"] and not ROOTLIKE_NAME in row["rgft_name"] and not ROOTLIKE_NAME in row["prune_name"]:
			features_restree_dict = calc_leaves_features(tree, "res", rgft_node_name=row["rgft_name"])  # res tree before optimization!
			df_prune = index_shared_features(df_prune, ind, row["prune_name"], "prune",  features_prune_dicts_dict)
			df_rgft = index_shared_features(df_rgft, ind, row["rgft_name"], "rgft", features_rgft_dicts_dict)
			df_rgft = index_additional_rgft_features(df_rgft, ind, row["prune_name"], row["rgft_name"], features_restree_dict, features_prune_dicts_dict)   # also prune dict because for 2 features i didn't want to comp dict within each rgft iteration (needed to compute on the starting tree)
			
			df_rgft.ix[ind, FEATURES["res_bl"]] = features_restree_dict['res_bl']
			df_rgft.ix[ind, FEATURES["res_tbl"]] = features_restree_dict['res_tbl']


	df_prune = df_prune[(df_prune["prune_name"] != ROOTLIKE_NAME) & (df_prune["rgft_name"] != ROOTLIKE_NAME)].dropna()
	df_rgft = df_rgft[(df_rgft["prune_name"] != ROOTLIKE_NAME) & (df_rgft["rgft_name"] != ROOTLIKE_NAME)].dropna()
	df_prune.to_csv(outpath_prune)  # runover existing one (with lls only) to fill in all features
	df_rgft.to_csv(outpath_rgft)    # runover existing one (with lls only) to fill in all features
	
	return

################################################################################################
######################### end of 'organizing features in csv' section ##########################
################################################################################################



if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='perform all SPR moves')
	parser.add_argument('--dataset_path', '-ds', required=True)
	args = parser.parse_args()

	dataset_path = args.dataset_path
	step_number = "1"
	outpath_prune = SUMMARY_PER_DS.format(dataset_path, "prune", "br", step_number)
	outpath_rgft = SUMMARY_PER_DS.format(dataset_path, "rgft", "br", step_number)
	
	collect_features(dataset_path, step_number, outpath_prune, outpath_rgft)