import sys

sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")

from defs import *
from utils.msa_functions import *
from data_processing.traverse_data_dirs import traverse_data_dirs
from summary.collect_SPR_features import *
from subprocess import Popen, PIPE, STDOUT
from parsing.parse_raxml_NG import parse_raxmlNG_content
from execute_programs.RAxML_NG import extract_model_params
from utils.create_job_file import get_job_qsub_command
import csv
from execute_programs import Phyml

ML_SOFTWARE_STARTING_TREE = 'phyml'     # could be phyml | RAxML_NG
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
	#if not nname:   # These 2 lines can be removed
	#	nname = "Nnew_p"

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
	# todo: verify what is the branch length in the pruning location (should be vp-v0 + vp-v1)

	return t_curr


def add_internal_names(tree_file, tree_file_cp_no_internal, t_orig):
	shutil.copy(tree_file, tree_file_cp_no_internal)
	for i, node in enumerate(t_orig.traverse()):
		if not node.is_leaf():
			node.name = "N{}".format(i)
	t_orig.write(format=3, outfile=tree_file)   # runover the orig file with no internal nodes names


def get_tree(ds_path, msa_file, rewrite_phylip, software=ML_SOFTWARE_STARTING_TREE):
	suf = "bionj" if not RANDOM_TREE_DIRNAME in ds_path else "br"
	tree_file = ds_path + PHYML_TREE_FILENAME.format(suf) if software == 'phyml' else ds_path + RAXML_TREE_FILENAME    # if software=='RAxML_NG'
	if rewrite_phylip:
		rewrite_in_phylip(msa_file)     # for one-time use on new ds

	#tree_file_cp_no_internal = ds_path + PHYML_TREE_FILENAME.format(suf + "_no_internal") if software == 'phyml' else ds_path + RAXML_TREE_FILENAME + "_no_internal"
	#if not os.path.exists(tree_file_cp_no_internal):
	#	t_orig = PhyloTree(newick=tree_file, alignment=msa_file, alg_format="iphylip", format=1)
	#	add_internal_names(tree_file, tree_file_cp_no_internal, t_orig)
	#else:
		#t_orig = PhyloTree(newick=tree_file, alignment=msa_file, alg_format="iphylip", format=3)
	t_orig = PhyloTree(newick=tree_file, alignment=msa_file, alg_format="iphylip", format=1)

	return t_orig


def call_phyml_storage(tree_dirpath, file_name, msa_file, runover, job_priority, cpmsa=False):
	opt_mode = 'br'  #if software == 'phyml' else 'fixed_subs'  # if software=='RAxML_NG'
	tree_path = tree_dirpath + file_name + ".txt"
	job_name = "phyml_" + "_".join([re.search("{}/*(.+?)/".format(DATA_PATH), tree_dirpath).group(1), tree_dirpath.split(SEP)[-3], file_name, opt_mode])

	cmd = "python " + CODE_PATH + "execute_programs/phyml.py " + "-f " + msa_file \
		  + " -br " + opt_mode + " -t " + tree_path
	if runover:
		cmd += " -r "
	if cpmsa:
		cmd += " -cp "

	os.system(cmd)


def call_raxml_mem(tree_str, msa_tmpfile, rates, pinv, alpha, freq):
	model_line_params = 'GTR{rates}+I{pinv}+G{alpha}+F{freq}'.format(rates="{{{0}}}".format("/".join(rates)),
									 pinv="{{{0}}}".format(pinv), alpha="{{{0}}}".format(alpha),
									 freq="{{{0}}}".format("/".join(freq)))

	# create tree file in memory and not in the storage:
	tree_rampath = "/dev/shm/" + msa_tmpfile.split(SEP)[-1] + "tree"  # the var is the str: tmp{dir_suffix}

	try:
		with open(tree_rampath, "w") as fpw:
			fpw.write(tree_str)

		p = Popen([RAXML_NG_SCRIPT, '--evaluate', '--msa', msa_tmpfile,'--threads', '1', '--opt-branches', 'on', '--opt-model', 'off', '--model', model_line_params, '--nofiles', '--tree', tree_rampath], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
		raxml_stdout = p.communicate()[0]
		raxml_output = raxml_stdout.decode()
		#print("\n"+raxml_output+"\n")

		res_dict = parse_raxmlNG_content(raxml_output)
		ll = res_dict['ll']
		rtime = res_dict['time']

	except Exception as e:
		print(msa_tmpfile.split(SEP)[-1][3:])
		print(e)
		exit()
	finally:
		os.remove(tree_rampath)

	return ll, rtime



def create_SPR_job(dataset_path, step_number, tree_type, rewrite_phy, runover):
	print("**************************************\n", dataset_path)
	job_name = "SPR_for_ds.sh"
	cmd = "python " + CODE_PATH + "execute_programs/SPR_move.py -ds " + dataset_path + " -st " + str(step_number) + " -ttype " + tree_type

	if runover:
		cmd += " -r "
	if rewrite_phy:
		cmd += " -phy "

	qsub_cmd = get_job_qsub_command(job_name=job_name,
									command=cmd,
									error_files_path=dataset_path + "error_files/")
	os.system(qsub_cmd)




def all_SPR(ds_path, outpath, tree=None, rewrite_phylip=False):
	orig_msa_file = ds_path + MSA_PHYLIP_FILENAME
	suf = "bionj" if not RANDOM_TREE_DIRNAME in ds_path else "br"
	stats_filepath = ds_path + PHYML_STATS_FILENAME.format(suf) if ML_SOFTWARE_STARTING_TREE == 'phyml' else ds_path + RAXML_STATS_FILENAME

	if 'ml_minus1' in ds_path:
		orig_msa_file = "/groups/itay_mayrose/danaazouri/PhyAI/DBset2/data/training_datasets/exampleSml/masked_species_real_msa.phy"
		stats_filepath = "/groups/itay_mayrose/danaazouri/PhyAI/DBset2/data/training_datasets/exampleSml/masked_species_real_msa.phy_phyml_stats_bionj.txt"

	t_orig = get_tree(ds_path, orig_msa_file, rewrite_phylip) if not tree else PhyloTree(newick=tree, alignment=orig_msa_file, alg_format="iphylip", format=1)
	if 'ml_minus1' in ds_path:
		t_orig.unroot()
		t_orig.get_tree_root().name = ROOTLIKE_NAME+"_2"
	else:
		t_orig.get_tree_root().name = ROOTLIKE_NAME if not tree else ROOTLIKE_NAME + "_2"


	st = "1" if not tree else "2"
	OUTPUT_TREES_FILE = TREES_PER_DS.format(ds_path, st)
	with open(OUTPUT_TREES_FILE, "w", newline='') as fpw:
		csvwriter = csv.writer(fpw)
		csvwriter.writerow(['', 'prune_name', 'rgft_name', 'newick'])

	# first, copy msa file to memory and save it:
	msa_rampath = "/dev/shm/tmp" + ds_path.split(SEP)[-2] #  to be on the safe side (even though other processes shouldn't be able to access it)

	with open(orig_msa_file) as fpr:
		msa_str = fpr.read()

	try:
		with open(msa_rampath, "w") as fpw:
			fpw.write(msa_str)  # don't write the msa string to a variable (or write and release it)
		msa_str = ''

		params_dict = (parse_phyml_stats_output(None, stats_filepath)) if ML_SOFTWARE_STARTING_TREE == 'phyml' else parse_raxmlNG_output(stats_filepath)
		freq, rates, pinv, alpha = [params_dict["fA"], params_dict["fC"], params_dict["fG"], params_dict["fT"]], [params_dict["subAC"], params_dict["subAG"], params_dict["subAT"], params_dict["subCG"],params_dict["subCT"], params_dict["subGT"]], params_dict["pInv"], params_dict["gamma"]
		df = pd.DataFrame()
		for i, prune_node in enumerate(t_orig.iter_descendants("levelorder")):
			prune_name = prune_node.name
			nname, subtree1, subtree2 = prune_branch(t_orig, prune_name) # subtree1 is the pruned subtree. subtree2 is the remaining subtree
			with open(OUTPUT_TREES_FILE, "a", newline='') as fpa:
				csvwriter = csv.writer(fpa)
				csvwriter.writerow([str(i)+",0", prune_name, SUBTREE1, subtree1.write(format=1)])
				csvwriter.writerow([str(i)+",1", prune_name, SUBTREE2, subtree2.write(format=1)])

			for j, rgft_node in enumerate(subtree2.iter_descendants("levelorder")):
				ind = str(i) + "," + str(j)
				rgft_name = rgft_node.name
				if nname == rgft_name: # if the rgrft node is the one that was pruned
					continue
				rearr_tree_str = regraft_branch(subtree2, rgft_node, subtree1, rgft_name, nname).write(format=1)

				### save tree to file by using "append"
				with open(OUTPUT_TREES_FILE, "a", newline='') as fpa:
					csvwriter = csv.writer(fpa)
					csvwriter.writerow([ind, prune_name, rgft_name, rearr_tree_str])

				#ll_rearr, rtime = call_raxml_mem(rearr_tree_str, msa_rampath, rates, pinv, alpha, freq)    # todo: uncomment ! temp for large ds

				df.loc[ind, "prune_name"], df.loc[ind, "rgft_name"] = prune_name, rgft_name
				#df.loc[ind, "time"] = rtime			# todo: uncomment ! temp for large ds
				#df.loc[ind, "ll"] = ll_rearr			# todo: uncomment ! temp for large ds


		df["orig_ds_ll"] = float(params_dict["ll"])
		df.to_csv(outpath.format("prune"))
		df.to_csv(outpath.format("rgft"))

	except Exception as e:
		print('could not complete the all_SPR function on dataset:', dataset_path, '\nError message:')
		print(e)
		exit()
	finally:
		os.remove(msa_rampath)
	return



# -st 1 -istart 0 -nlines 1   ## -phy {if first time running anything on the dataset}
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='perform all SPR moves')
	parser.add_argument('--dataset_path', '-ds', default=None)
	parser.add_argument('--runover', '-r', default=False, action='store_true')
	parser.add_argument('--rewrite_in_phylip', '-phy', default=False, action='store_true')
	parser.add_argument('--tree_type', '-ttype', default='bionj') ### does mean NOTHING when runing on RAXML_NG mode    ## could be bionj or random
	parser.add_argument('--index_to_start_run', '-istart', default=False)
	parser.add_argument('--nline_to_run', '-nlines', default=False)		 # number of lines from the dataset (int): n^2
	parser.add_argument('--step_number', '-st', required=True)			 # counting from 1
	args = parser.parse_args()

	dataset_path = args.dataset_path
	if dataset_path:
		dataset_path = dataset_path if args.tree_type == 'bionj' else dataset_path + RANDOM_TREE_DIRNAME # if == 'random
		outpath = SUMMARY_PER_DS.format(dataset_path, "{}", "br", args.step_number)
		if args.runover:
			if args.step_number == "1":
				res = all_SPR(dataset_path, outpath, tree=None, rewrite_phylip=args.rewrite_in_phylip)
			else:   # run next step on previous' best tree
				prev_step = str(int(args.step_number)-1)
				dfr = pd.read_csv(TREES_PER_DS.format(dataset_path, prev_step), index_col=0)
				df_sum = pd.read_csv(SUMMARY_PER_DS.format(dataset_path, "prune", "br", prev_step)).set_index('Unnamed: 0')
				best_tree_id = df_sum["ll"].astype(float).idxmax()
				tree_str = dfr.loc[best_tree_id, "newick"]

				with open('/groups/itay_mayrose/danaazouri/PhyAI/DBset2/data/training_datasets/exampleSphaero/masked_species_real_msa.phy_phyml_tree_bionj_st2.txt', 'r') as fp:
					tree_str = fp.read().strip()
				res = all_SPR(dataset_path, outpath, tree=tree_str, rewrite_phylip=args.rewrite_in_phylip)

		collect_features(dataset_path, args.step_number, outpath.format("prune"), outpath.format("rgft"), args.tree_type)
	else:
		csv_path = SUMMARY_FILES_DIR + CHOSEN_DATASETS_FILENAME
		res = traverse_data_dirs(create_SPR_job, csv_path, (args.index_to_start_run, args.nline_to_run), args.step_number, args.tree_type, args.rewrite_in_phylip, args.runover)