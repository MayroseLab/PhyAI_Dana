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
############### begining of 'parsing rearrangements and PhyML outputs' section #################
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


def get_newick_tree(tree):
	"""
	:param tree: newick tree string or txt file containing one tree
	:return:	tree: a string of the tree in ete3.Tree format
	"""
	if os.path.exists(tree):
		with open(tree, 'r') as tree_fpr:
			tree = tree_fpr.read().strip()
	return tree


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
		print(treepath_no_internal)

	with open(treepath_no_internal, 'r') as fp:
		tree_str = fp.read()
		if re.search("\)N", tree_str):
			return
	with open(treepath_no_internal, 'w') as fp:
		fp.write(new_str)
		
		
def parse_phyml_stats_output(msa_filepath, stats_filepath):
	"""
	:return: dictionary with the attributes - string typed. if parameter was not estimated, empty string
	"""
	res_dict = dict.fromkeys(["ntaxa", "nchars", "ll",
	                          "fA", "fC", "fG", "fT",
	                          "subAC", "subAG", "subAT", "subCG", "subCT", "subGT",
	                          "pInv", "gamma",
	                          "path"], "")
	
	if msa_filepath:
		res_dict['ntaxa'], res_dict['nchars'] = (str(x) for x in get_msa_properties(get_msa_from_file(msa_filepath)))
	
	res_dict["path"] = stats_filepath
	try:
		with open(stats_filepath) as fpr:
			content = fpr.read()
		
		# likelihood
		res_dict["ll"] = re.search("Log-likelihood:\s+(.*)", content).group(1).strip()
		
		# gamma (alpha parameter) and proportion of invariant sites
		gamma_regex = re.search("Gamma shape parameter:\s+(.*)", content)
		pinv_regex = re.search("Proportion of invariant:\s+(.*)", content)
		if gamma_regex:
			res_dict['gamma'] = gamma_regex.group(1).strip()
		if pinv_regex:
			res_dict['pInv'] = pinv_regex.group(1).strip()
		
		# Nucleotides frequencies
		for nuc in "ACGT":
			nuc_freq = re.search("  - f\(" + nuc + "\)\= (.*)", content).group(1).strip()
			res_dict["f" + nuc] = nuc_freq
		
		# substitution frequencies
		for nuc1 in "ACGT":
			for nuc2 in "ACGT":
				if nuc1 < nuc2:
					nuc_freq = re.search(nuc1 + " <-> " + nuc2 + "(.*)", content).group(1).strip()
					res_dict["sub" + nuc1 + nuc2] = nuc_freq
	except:
		print("Error with:", res_dict["path"], res_dict["ntaxa"], res_dict["nchars"])
		return
	return res_dict


def return_ll(tree_dirpath, msa_file, filename, br_mode):
	stats_filepath = SEP.join([tree_dirpath, "{}_phyml_{}_{}.txt".format(filename, "stats", br_mode)])
	try:
		res_dict = parse_phyml_stats_output(msa_file, stats_filepath)
		ll_rearr = float(res_dict["ll"])
	except:
		ll_rearr = None
		pass
		#print("does not exist or empty")

	return ll_rearr


def parse_neighbors_dirs(ds_path, outpath_prune, outpath_rgft, step_number, cp_internal=False, tree_type="bionj"):
	'''
	this function is used only when re-running SPR. otherwise, I need to parse the 'newick.csv'
	'''
	print("**** ", ds_path)

	msa_file = ds_path + MSA_PHYLIP_FILENAME
	all_trees = ds_path + REARRANGEMENTS_NAME + "s/"
	outpath_trees = TREES_PER_DS.format(ds_path, step_number)

	res_dict_orig_tree = parse_phyml_stats_output(msa_file, ds_path + PHYML_STATS_FILENAME.format(tree_type))
	ll_orig_tree = float(res_dict_orig_tree["ll"])

	df = pd.DataFrame(index=np.arange(0))
	df2 = pd.DataFrame(index=np.arange(0))
	for i, prune_name in enumerate(os.listdir(all_trees)):
		prune_dirpath = SEP.join([all_trees, prune_name, ""])

		for j, rgft_name in enumerate(os.listdir(prune_dirpath)):
			ind = str(i) + "," +str(j)
			tree_dirpath = SEP.join([all_trees, prune_name, rgft_name, ""])
			treename = SUBTREE1 if j == 0 else SUBTREE2 if j == 1 else REARRANGEMENTS_NAME
			tree_path = SEP.join([tree_dirpath, "{}.txt".format(treename)])
			# save all rearrangements on one file per ds (before deleting all)
			df2.ix[ind, "prune_name"], df2.ix[ind, "rgft_name"] = prune_name, rgft_name
			df2.ix[ind, "newick"] = get_newick_tree(tree_path)

			if not "subtree" in rgft_name:  # subtrees are dealt separately ~10 lines above
				if cp_internal:
					treepath_with_internal = SEP.join([tree_dirpath, REARRANGEMENTS_NAME + ".txt"])
					rearr_tree_path = SEP.join([tree_dirpath, "{}_phyml_{}_{}.txt".format(MSA_PHYLIP_FILENAME, "tree", "br")])
					cp_internal_names(rearr_tree_path, treepath_with_internal)
				ll_rearr = return_ll(tree_dirpath, msa_file, MSA_PHYLIP_FILENAME, "br")

				df.ix[ind, "prune_name"], df.ix[ind, "rgft_name"] = prune_name, rgft_name
				df.ix[ind, "orig_ds_ll"], df.ix[ind, "ll"] = ll_orig_tree, ll_rearr

	df.to_csv(outpath_prune)
	df.to_csv(outpath_rgft)
	df2.to_csv(outpath_trees)

	# to reduce inodes number- delete subdirs after copying important content to 2 csvs:
	shutil.rmtree(all_trees, ignore_errors=True)

	return

################################################################################################
################## end of 'parsing rearrangements and PhyML outputs' section ###################
################################################################################################
################################################################################################
######################### begining of 'generate SPR neigbors' section ##########################
################################################################################################

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


def regraft_branch(t_cp_p, rgft_node, prune_node_cp, rgft_name, nname):
	'''
	get a tree with the 2 concatenated subtrees
	'''

	new_branch_length = rgft_node.dist /2
	t_temp = PhyloTree()  			   # for concatenation of both subtrees ahead, to avoid polytomy
	t_temp.add_child(prune_node_cp)
	t_curr = t_cp_p.copy()
	rgft_node_cp = t_curr & rgft_name  # locate the node in the copied subtree

	rgft_loc = rgft_node_cp.up
	rgft_node_cp.detach()
	t_temp.add_child(rgft_node_cp, dist=new_branch_length)
	t_temp.name = nname
	rgft_loc.add_child(t_temp, dist=new_branch_length)  # regrafting

	return t_curr


def add_internal_names(tree_file, tree_file_cp_no_internal, t_orig):
	shutil.copy(tree_file, tree_file_cp_no_internal)
	for i, node in enumerate(t_orig.traverse()):
		if not node.is_leaf():
			node.name = "N{}".format(i)
	t_orig.write(format=3, outfile=tree_file)   # runover the orig file with no internal nodes names


def get_tree(ds_path, msa_file):
	tree_file = ds_path + PHYML_TREE_FILENAME.format("bionj")
	tree_file_cp_no_internal = ds_path + PHYML_TREE_FILENAME.format("bionj_no_internal")
	if not os.path.exists(tree_file_cp_no_internal):
		t_orig = PhyloTree(newick=tree_file, alignment=msa_file, alg_format="iphylip", format=1)
		add_internal_names(tree_file, tree_file_cp_no_internal, t_orig)
	else:
		t_orig = PhyloTree(newick=tree_file, alignment=msa_file, alg_format="iphylip", format=3)

	return t_orig


def save_rearr_file(trees_dirpath, rearrtree, filename, runover=False):
	if not os.path.exists(trees_dirpath):
		os.makedirs(trees_dirpath)
	tree_path = trees_dirpath + filename + ".txt"
	if runover or not os.path.exists(tree_path):
		rearrtree.write(format=1, outfile=tree_path)

	return tree_path


def call_phyml(tree_dirpath, file_name, msa_file, runover, job_priority, br_mode, cpmsa=False):
	tree_path = tree_dirpath + file_name + ".txt"
	job_name = "phyml_" + "_".join([re.search("{}/*(.+?)/".format(DATA_PATH), tree_dirpath).group(1), tree_dirpath.split(SEP)[-3], file_name, br_mode])
	cmd = "python " + PHYML_SCRIPT_PATH + " -f " + msa_file \
		  + " -br " + br_mode + " -t " + tree_path
	if runover:
		cmd += " -r "
	if cpmsa:
		cmd += " -cp "

	os.system(cmd)


def all_SPR(ds_path, tree=None, runover=False, job_priority=-1):
	orig_msa_file = ds_path + MSA_PHYLIP_FILENAME
	t_orig = get_tree(ds_path, orig_msa_file) if not tree else PhyloTree(newick=tree, alignment=orig_msa_file, alg_format="iphylip", format=1)
	t_orig.get_tree_root().name = ROOTLIKE_NAME if not tree else ROOTLIKE_NAME+"_2"
	for i, prune_node in enumerate(t_orig.iter_descendants("levelorder")):
		prune_name = prune_node.nam
		nname, subtree1, subtree2 = prune_branch(t_orig, prune_name) # subtree1 is the pruned subtree. subtree2 is the remaining subtree
		subtrees_dirpath = SEP.join([ds_path, REARRANGEMENTS_NAME+"s", prune_name, "{}", ""])

		for j, rgft_node in enumerate(subtree2.iter_descendants("levelorder")):
			rgft_name = rgft_node.name
			if nname == rgft_name: # if the rgrft node is the one that was pruned
				continue

			full_tree_dirpath = subtrees_dirpath.format(rgft_name)
			if runover or not os.path.exists(full_tree_dirpath + REARRANGEMENTS_NAME + ".txt"):
				full_tree = regraft_branch(subtree2, rgft_node, subtree1, rgft_name, nname)
				save_rearr_file(full_tree_dirpath, full_tree, filename=REARRANGEMENTS_NAME, runover=runover)
				call_phyml(full_tree_dirpath, REARRANGEMENTS_NAME, orig_msa_file, runover, job_priority, "br", cpmsa=True)
	return

################################################################################################
########################### end of 'generate SPR neigbors' section #############################
################################################################################################

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='perform all SPR moves')
	parser.add_argument('--dataset_path', '-ds', required=True)
	args = parser.parse_args()
	
	dataset_path = args.dataset_path
	outpath_prune = SUMMARY_PER_DS.format(dataset_path, "prune", "br", "")
	outpath_rgft = SUMMARY_PER_DS.format(dataset_path, "rgft", "br", step_number="1")
	
	if not os.path.exists(outpath_rgft) or not os.path.exists(outpath_prune):
		all_SPR(dataset_path, tree=None, runover=False)
		parse_neighbors_dirs(dataset_path, outpath_prune, outpath_rgft, step_number="1", cp_internal=True, tree_type='bionj')