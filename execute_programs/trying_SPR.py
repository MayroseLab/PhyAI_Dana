import sys

sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")

from defs import *
from ete3 import Tree, PhyloTree
from execute_programs.Phyml import run_phyml
from parsing.parse_phyml import parse_phyml_stats_output
from utils.msa_functions import *
from data_processing.traverse_data_dirs import traverse_data_dirs
from summary.collect_SPR_features import *
from random import randint
from subprocess import Popen, PIPE, STDOUT, run
from parsing.parse_raxml_NG import parse_raxmlNG_content
from execute_programs.RAxML_NG import extract_model_params


ML_SOFTWARE = 'RAxML_NG'     # could be phyml | RAxML_NG
RAXML_NG_SCRIPT = "raxml-ng"




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
	# todo: verify what is the branch length in the pruning location (should be vp-v0 + vp-v1)

	return t_curr


def add_internal_names(tree_file, tree_file_cp_no_internal, t_orig):
	shutil.copy(tree_file, tree_file_cp_no_internal)
	for i, node in enumerate(t_orig.traverse()):
		if not node.is_leaf():
			node.name = "N{}".format(i)
	t_orig.write(format=3, outfile=tree_file)   # runover the orig file with no internal nodes names


def get_tree(ds_path, msa_file, rewrite_phylip):
	#suf = "bionj" if not RANDOM_TREE_DIRNAME in ds_path else "br"
	suf = "bionj" if not RANDOM_TREE_DIRNAME in ds_path else "br"
	tree_file = ds_path + PHYML_TREE_FILENAME.format(suf) if ML_SOFTWARE == 'phyml' else ds_path + RAXML_TREE_FILENAME    # if software=='RAxML_NG'
	if rewrite_phylip:
		rewrite_in_phylip(msa_file)     # for one-time use on new ds
	tree_file_cp_no_internal = ds_path + PHYML_TREE_FILENAME.format(suf + "_no_internal") if ML_SOFTWARE == 'phyml' else ds_path + RAXML_TREE_FILENAME + "_no_internal"
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


def call_ml_software(tree_dirpath, file_name, msa_file, runover, job_priority, cpmsa=False, software=ML_SOFTWARE):
	opt_mode = 'br' if software == 'phyml' else 'fixed_subs'  # if software=='RAxML_NG'
	tree_path = tree_dirpath + file_name + ".txt"
	job_name = "{}_".format(software) + "_".join([re.search("{}/*(.+?)/".format(DATA_PATH), tree_dirpath).group(1), tree_dirpath.split(SEP)[-3], file_name, opt_mode])

	cmd = "python " + CODE_PATH + "execute_programs/{}.py ".format(software) + "-f " + msa_file \
		  + " -br " + opt_mode + " -t " + tree_path
	if runover:
		cmd += " -r "
	if cpmsa:
		cmd += " -cp "

	os.system(cmd)

def call_raxml_mem(tree_str, msa_file, rates, pinv, alpha, freq):
	model_line_params = 'GTR{rates}+I{pinv}+G{alpha}+F{freq}'.format(rates="{{{0}}}".format("/".join(rates)),
									 pinv="{{{0}}}".format(pinv), alpha="{{{0}}}".format(alpha),
									 freq="{{{0}}}".format("/".join(freq)))

	# now prepare the raxml command, "--tree" should be last because the function behaves as the pipe
	p = run([RAXML_NG_SCRIPT, '--evaluate', '--msa', msa_file,'--threads', '1', '--opt-branches', 'on', '--opt-model', 'off', '--model', model_line_params, '--nofiles', '--tree'],
			stdout=PIPE, input=tree_str, encoding='ascii')

	res_dict = parse_raxmlNG_content(p.stdout)
	ll = res_dict['ll']

	return ll


def remove_redundant_nodes(tree, ntaxa, node_name = None):
	if len(tree.get_descendants()) > ntaxa * 2 - 3:
		tree.get_tree_root().children[0].delete(preserve_branch_length=True)
	#if node_name and len(tree.get_descendants()) > ntaxa * 2 - 3:
	#	tree.search_nodes(name=node_name)[0].up.delete(preserve_branch_length=True)

	return tree


def create_SPR_job(dataset_path, step_number, tree_type, rewrite_phy, runover):
	print("**************************************\n", dataset_path)
	job_name = "SPR_for_ds"
	cmd = "python " + CODE_PATH + "execute_programs/SPR_move.py -ds " + dataset_path + " -st " + str(step_number) + " -ttype " + tree_type
	if runover:
		cmd += " -r "
	if rewrite_phy:
		cmd += " -phy "
	create_job_file.main(command=cmd, dirpath=dataset_path, sh_file=job_name + ".sh", multiply_jobs=False,
						 priority=-1, job_name=job_name)





def all_SPR(ds_path, tree=None, rewrite_phylip=False, runover=False, job_priority=-1):
	orig_msa_file = ds_path + MSA_PHYLIP_FILENAME
	t_orig = get_tree(ds_path, orig_msa_file, rewrite_phylip) if not tree else PhyloTree(newick=tree, alignment=orig_msa_file, alg_format="iphylip", format=1)
	t_orig.get_tree_root().name = ROOTLIKE_NAME if not tree else ROOTLIKE_NAME+"_2"
	OUTPUT_TREES_FILE = "trees.txt"
	with open(OUTPUT_TREES_FILE, "w") as fpw:
		fpw.write(",prune_name,rgft_name,newick")

	# first, copy msa file to memory and save it:
	msa_rampath = "/dev/shm/tmp" + str(randint(10 ** 6,
											   10 ** 7 - 1))  # random number to be on the safe side (even though other processes shouldn't be able to access it)
	with open(msa_rampath, "w") as fpw:
		with open(orig_msa_file) as fpr:
			fpw.write(fpr.read())  # don't write the msa string to a variable (or write and release it)

	rates, pinv, alpha, freq = extract_model_params(msa_rampath, ds_path, 'raxml')
	df = pd.DataFrame()
	for i, prune_node in enumerate(t_orig.iter_descendants("levelorder")):
		prune_name = prune_node.name
		nname, subtree1, subtree2 = prune_branch(t_orig, prune_name) # subtree1 is the pruned subtree. subtree2 is the remaining subtree
		subtrees_dirpath = SEP.join([ds_path, REARRANGEMENTS_NAME+"s", prune_name, "{}", ""])
		#call_phyml_ll(subtrees_dirpath, subtree1, subtree2, seqs_dict, runover=runover) # cal phyml foreach subtree + truncated msa
		#save_rearr_file(subtrees_dirpath.format(SUBTREE1), subtree1, filename=SUBTREE1,runover=runover)  # will print it out WITHOUT the root NAME !
		#save_rearr_file(subtrees_dirpath.format(SUBTREE2), subtree2, filename=SUBTREE2,runover=runover)  # will print it out WITHOUT the root NAME !
		# todo: continue
		#with open(OUTPUT_TREES_FILE, "a") as fpa:
		#	fpa.write(", ".join([ind, prune_name, rgft_name, rearr_tree_str]))


		for j, rgft_node in enumerate(subtree2.iter_descendants("levelorder")):
			ind = str(i) + "," + str(j)
			rgft_name = rgft_node.name
			if nname == rgft_name: # if the rgrft node is the one that was pruned
				continue

			rearr_tree_str = regraft_branch(subtree2, rgft_node, subtree1, rgft_name, nname).write(format=1)
			### save tree to file by using "append"
			with open(OUTPUT_TREES_FILE, "a") as fpa:
				fpa.write(", ".join([ind, prune_name, rgft_name, rearr_tree_str]))


			ll_rearr = call_raxml_mem(rearr_tree_str, msa_rampath, rates, pinv, alpha, freq)
			print(ll_rearr)
			exit()
			df.loc[ind, "prune_name"], df.loc[ind, "rgft_name"] = prune_name, rgft_name
			df.loc[ind, "ll"] = ll_rearr
		#exit()
	df["orig_ds_ll"] = float(parse_raxmlNG_output(ds_path + RAXML_STATS_FILENAME)['ll'])
	df.to_csv(SUMMARY_PER_DS.format(dataset_path, "lls", 'br', '1'))   # todo: update
	return



# -phy -cp -istart 0 -nlines 1 -st 1
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='perform all SPR moves')
	parser.add_argument('--dataset_path', '-ds', default=None)
	parser.add_argument('--runover', '-r', default=False, action='store_true')
	parser.add_argument('--rewrite_in_phylip', '-phy', default=False, action='store_true')
	parser.add_argument('--cp_internal', '-cp', default=False, action='store_true')
	parser.add_argument('--tree_type', '-ttype', default='bionj')  # could be bionj or random
	parser.add_argument('--index_to_start_run', '-istart', default=False)
	parser.add_argument('--nline_to_run', '-nlines', default=False)		 # number of lines from the dataset (int): n^2
	parser.add_argument('--step_number', '-st', required=True)			 # counting from 1
	args = parser.parse_args()

	dataset_path = args.dataset_path
	if dataset_path:
		dataset_path = dataset_path if args.tree_type == 'bionj' else dataset_path + RANDOM_TREE_DIRNAME # if == 'random
		outpath_prune = SUMMARY_PER_DS.format(dataset_path, "prune", "br", args.step_number)
		outpath_rgft = SUMMARY_PER_DS.format(dataset_path, "rgft", "br", args.step_number)
		# todo: uncomment the folowwing condition after running over the failed runings
		#if not os.path.exists(outpath_rgft) or not os.path.exists(outpath_prune):
		if args.step_number == "1":
			res = all_SPR(dataset_path, tree=None, rewrite_phylip=args.rewrite_in_phylip, runover=args.runover, job_priority=-1)
		else:   # run next step on previous' best tree
			prev_step = str(int(args.step_number)-1)
			dfr = pd.read_csv(TREES_PER_DS.format(dataset_path, prev_step), index_col=0)
			df_sum = pd.read_csv(SUMMARY_PER_DS.format(dataset_path, "prune", "br", prev_step)).set_index('Unnamed: 0')
			best_tree_id = df_sum["ll"].astype(float).idxmax()
			tree_str = dfr.loc[best_tree_id, "newick"]

			res = all_SPR(dataset_path, tree=tree_str, rewrite_phylip=args.rewrite_in_phylip, runover=args.runover, job_priority=-1)

		# after preforming all steps, parse results and delete all files
		res = parse_neighbors_dirs(dataset_path, outpath_prune, outpath_rgft, args.step_number, args.cp_internal, tree_type=args.tree_type)
		collect_features(dataset_path, args.step_number, outpath_prune, outpath_rgft, args.tree_type)
	else:
		csv_path = SUMMARY_FILES_DIR + CHOSEN_DATASETS_FILENAME
		res = traverse_data_dirs(create_SPR_job, csv_path, (args.index_to_start_run, args.nline_to_run), args.step_number, args.tree_type, args.rewrite_in_phylip, args.runover)