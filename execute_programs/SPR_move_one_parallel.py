import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")

from defs import *
from ete3 import Tree, PhyloTree
from execute_programs.Phyml import run_phyml
from parsing.parse_phyml import parse_phyml_stats_output
from utils.msa_functions import *
from data_processing.traverse_data_dirs import traverse_data_dirs
from summary.collect_SPR_features import *



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
	# print("### rgft node: ", rgft_name)

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





def get_tree(ds_path, msa_file, rewrite_phylip):
	suf = "bionj" if not RANDOM_TREE_DIRNAME in ds_path else "br"
	tree_file = ds_path + PHYML_TREE_FILENAME.format(suf)
	if rewrite_phylip:
		rewrite_in_phylip(msa_file)     # for one-time use on new ds
	tree_file_cp_no_internal = ds_path + PHYML_TREE_FILENAME.format(suf + "_no_internal")
	#print("\n################ tree file:" + tree_file + " ################\n")
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
	cmd = "python " + CODE_PATH + "execute_programs/Phyml.py " + "-f " + msa_file \
		  + " -br " + br_mode + " -t " + tree_path
	if runover:
		cmd += " -r "
	if cpmsa:
		cmd += " -cp "

	#os.system(cmd)
	create_job_file.main(command=cmd, dirpath=tree_dirpath, sh_file=job_name + ".sh", multiply_jobs=False,
	                     priority=-1, job_name=job_name)



def remove_redundant_nodes(tree, ntaxa, node_name = None):
	if len(tree.get_descendants()) > ntaxa * 2 - 3:
		tree.get_tree_root().children[0].delete(preserve_branch_length=True)
	#if node_name and len(tree.get_descendants()) > ntaxa * 2 - 3:
	#	tree.search_nodes(name=node_name)[0].up.delete(preserve_branch_length=True)

	return tree



def subtrees(subtrees_dirpath, subtree1, subtree2, runover=False, job_priority=-1):
	'''
	First truncate msa according to each subtree. Then call phyml with the respective tree+trunc_msa file
	'''
	for subtree_name in [SUBTREE1, SUBTREE2]:  									  #do same for both subtrees result from prunning
		tree = subtree1 if subtree_name == SUBTREE1 else subtree2
		subtree_dirpath = subtrees_dirpath.format(subtree_name)
		if not os.path.exists(subtree_dirpath):
			os.makedirs(subtree_dirpath)

		treepath = save_rearr_file(subtree_dirpath, tree, filename=subtree_name, runover=runover)  # will print it out WITHOUT the root NAME !



def create_SPR_job(dataset_path, step_number, tree_type, rewrite_phy, runover):
	print("**************************************\n", dataset_path)
	job_name = "SPR_for_ds"
	cmd = "python " + CODE_PATH + "execute_programs/SPR_move_one_parallel.py -ds " + dataset_path + " -st " + str(step_number) + " -ttype " + tree_type
	if runover:
		cmd += " -r "
	if rewrite_phy:
		cmd += " -phy "
	create_job_file.main(command=cmd, dirpath=dataset_path, sh_file=job_name + ".sh", multiply_jobs=False,
						 priority=-1, job_name=job_name)
	
	

def create_SPR_i_job(dataset_path, orig_msa_file, prune_name, tree_type, rewrite_phy, runover):
	job_name = "SPR_i"
	cmd = "python " + CODE_PATH + "execute_programs/SPR_move_one_parallel.py -ds " + dataset_path + " -pr " + prune_name + " -ttype " + tree_type + " -msa " + orig_msa_file
	if runover:
		cmd += " -r "
	if rewrite_phy:
		cmd += " -phy "
	create_job_file.main(command=cmd, dirpath=dataset_path, sh_file=job_name + ".sh", multiply_jobs=False,
						 priority=-1, job_name=job_name)




def all_SPR_i(ds_path, orig_msa_file, prune_name, runover=False, job_priority=-1):
	t_orig = get_tree(dataset_path, orig_msa_file, args.rewrite_in_phylip)
	t_orig.get_tree_root().name = ROOTLIKE_NAME
	
	nname, subtree1, subtree2 = prune_branch(t_orig, prune_name) # subtree1 is the pruned subtree. subtree2 is the remaining subtree
	subtrees_dirpath = SEP.join([ds_path, REARRANGEMENTS_NAME+"s", prune_name, "{}", ""])
	subtrees(subtrees_dirpath, subtree1, subtree2,runover=runover)

	for j, rgft_node in enumerate(subtree2.iter_descendants("levelorder")):
		rgft_name = rgft_node.name
		if nname == rgft_name: # if the rgrft node is the one that was pruned
			continue
		
		full_tree_dirpath = subtrees_dirpath.format(rgft_name)
		
		
		f = SEP.join([full_tree_dirpath, "{}_phyml_{}_{}.txt".format(MSA_PHYLIP_FILENAME, "stats", "br")])
		try:
			with open(f) as fpr:
				content = fpr.read()
				ll_rearr = float(re.search("Log-likelihood:\s+(.*)", content).group(1).strip())
		except:
			print(f)
			ll_rearr = False
		

		
	
		if runover or not os.path.exists(full_tree_dirpath + REARRANGEMENTS_NAME + ".txt") or not ll_rearr:
			full_tree = regraft_branch(subtree2, rgft_node, subtree1, rgft_name, nname)
			#remove_redundant_nodes(full_tree, ntaxa, rgft_name)
			save_rearr_file(full_tree_dirpath, full_tree, filename=REARRANGEMENTS_NAME, runover=runover)
			call_phyml(full_tree_dirpath, REARRANGEMENTS_NAME, orig_msa_file, runover, job_priority, "br", cpmsa=True)
		else:
			print(full_tree_dirpath + " already exists")
	#exit()

	return



# -phy -cp -istart 0 -nlines 1 -st 1
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='perform all SPR moves')
	parser.add_argument('--dataset_path', '-ds', default=None)
	parser.add_argument('--runover', '-r', default=False, action='store_true')
	parser.add_argument('--rewrite_in_phylip', '-phy', default=False, action='store_true')
	parser.add_argument('--cp_internal', '-cp', default=False, action='store_true')
	parser.add_argument('--tree_type', '-ttype', default='bionj')  # could be bionj or random
	parser.add_argument('--msa_path', '-msa', default=False)
	parser.add_argument('--prune_name', '-pr', default=False)
	args = parser.parse_args()
	
	dataset_path = args.dataset_path
	if dataset_path:
		dataset_path = dataset_path if args.tree_type == 'bionj' else dataset_path + RANDOM_TREE_DIRNAME # if == 'random
		outpath_prune = SUMMARY_PER_DS.format(dataset_path, "prune", "br", "1")
		outpath_rgft = SUMMARY_PER_DS.format(dataset_path, "rgft", "br", "1")
		if not os.path.exists(outpath_rgft) or not os.path.exists(outpath_prune):
			res = all_SPR_i(dataset_path, args.msa_path, args.prune_name, runover=args.runover, job_priority=-1)
	
	else:
		dataset_path = DATA_PATH + "exampleV2"
		orig_msa_file = dataset_path + MSA_PHYLIP_FILENAME
		t_orig = get_tree(dataset_path, orig_msa_file, args.rewrite_in_phylip)
		t_orig.get_tree_root().name = ROOTLIKE_NAME
		
		for i, prune_node in enumerate(t_orig.iter_descendants("levelorder")):
			if i <= 1:
				prune_name = prune_node.name
				create_SPR_i_job(dataset_path, orig_msa_file, prune_name, args.tree_type, args.rewrite_in_phylip, args.runover)