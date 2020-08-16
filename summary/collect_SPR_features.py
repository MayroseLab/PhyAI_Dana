import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")

from defs import *
from parsing.parse_phyml import parse_phyml_stats_output
from utils.tree_functions import *
from data_processing.traverse_data_dirs import traverse_data_dirs

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
from utils.tree_functions import *
from execute_programs.SPR_move import *

OPT_TYPE = "br"


def index_additional_rgft_features(df_rgft, ind, prune_name, rgft_name, features_restree_dict, features_dict_prune):
	df_rgft.ix[ind, FEATURES["top_dist"]] = features_dict_prune['top_dist'][prune_name][rgft_name]
	df_rgft.ix[ind, FEATURES["bl_dist"]] = features_dict_prune['bl_dist'][prune_name][rgft_name]
	# updated late (23/4)
	df_rgft.ix[ind, FEATURES["res_bl"]] = features_restree_dict['res_bl']
	df_rgft.ix[ind, FEATURES["res_tbl"]] = features_restree_dict['res_tbl']

	return df_rgft




def return_ll(tree_dirpath, br_mode):
	msa_filepath = tree_dirpath + MSA_PHYLIP_FILENAME
	stats_filepath = "{}_phyml_{}_{}.txt".format(msa_filepath, "stats", br_mode)
	try:
		res_dict = parse_phyml_stats_output(msa_filepath, stats_filepath)
		ll_rearr = float(res_dict["ll"])
	except:
		ll_rearr = None
		pass
		#print("does not exist or empty")

	return ll_rearr




def index_shared_features(dff, ind, edge, move_type, features_dicts_dict):
	d_ll = format(dff.loc[ind, "ll"] - dff.loc[ind, "orig_ds_ll"], '.20f')
	dff.ix[ind, LABEL.format(move_type)] = d_ll  # LABEL

	#*tbl of orig tree will be calculated via 'RandomForest_learning' script		#f1  (FeaturePruning)
	dff.ix[ind, FEATURES["bl"]] = features_dicts_dict["bl"][edge] 					#f2
	dff.ix[ind, FEATURES["longest"]] = features_dicts_dict["longest"]				#f3


	#index subtrees features
	for subtype in ["p", "r"]:
		dff.ix[ind, FEATURES["ntaxa_{}".format(subtype)]] = features_dicts_dict["ntaxa_{}".format(subtype)][edge]  		#f4,5
		#dff.ix[ind, FEATURES["pdist_{}".format(subtype)]] = features_dicts_dict["pdist_{}".format(subtype)][edge]    	#f6,7
		dff.ix[ind, FEATURES["tbl_{}".format(subtype)]] = features_dicts_dict["tbl_{}".format(subtype)][edge]  			#f8,9
		#dff.ix[ind, FEATURES["pars_{}".format(subtype)]] = features_dicts_dict["pars_{}".format(subtype)][edge] 		#f10,11
		dff.ix[ind, FEATURES["longest_{}".format(subtype)]] = features_dicts_dict["longest_{}".format(subtype)][edge]	#f12,13

	return dff



def collect_features(ds_path, step_number, outpath_prune, outpath_rgft, tree_type='bionj'):
	dfr = pd.read_csv(TREES_PER_DS.format(ds_path, step_number), index_col=0)
	orig_ds_msa_file = ds_path + MSA_PHYLIP_FILENAME
	df_prune = pd.read_csv(outpath_prune, index_col=0)
	df_rgft = pd.read_csv(outpath_rgft, index_col=0)

	suf = "bionj" if tree_type == 'bionj' else 'br'  # if tree_type="random"
	features_prune_dicts_dict = calc_leaves_features(ds_path + PHYML_TREE_FILENAME.format(suf),"prune")
	#features_prune_dicts_dict = calc_leaves_features(ds_path + PHYML_TREE_FILENAME.format(suf), orig_ds_msa_file, "prune")

	for i, row in dfr.iterrows():
		ind = row.name
		print(ind)
		tree = row["newick"]
		if row["rgft_name"] == SUBTREE2:	# namely the remaining subtree
			features_rgft_dicts_dict = calc_leaves_features(tree, "rgft")
			#features_rgft_dicts_dict = calc_leaves_features(tree, orig_ds_msa_file, "rgft")  # msa will be truncated within the function
			#pass
		if not "subtree" in row["rgft_name"] and not ROOTLIKE_NAME in row["rgft_name"] and not ROOTLIKE_NAME in row["prune_name"]:
			features_restree_dict = calc_leaves_features(tree, "res", rgft_node_name=row["rgft_name"])
			#features_restree_dict = calc_leaves_features(tree, orig_ds_msa_file, "res", rgft_node_name=row["rgft_name"])
			df_prune = index_shared_features(df_prune, ind, row["prune_name"], "prune",  features_prune_dicts_dict)
			df_rgft = index_shared_features(df_rgft, ind, row["rgft_name"], "rgft", features_rgft_dicts_dict)
			df_rgft = index_additional_rgft_features(df_rgft, ind, row["prune_name"], row["rgft_name"], features_restree_dict, features_prune_dicts_dict)   # also prune dict because for 2 features i didn't want to comp dict within each rgft iteration (needed to compute on the starting tree)
			
			df_rgft.ix[ind, FEATURES["res_bl"]] = features_restree_dict['res_bl']
			df_rgft.ix[ind, FEATURES["res_tbl"]] = features_restree_dict['res_tbl']


	df_prune = df_prune[(df_prune["prune_name"] != ROOTLIKE_NAME) & (df_prune["rgft_name"] != ROOTLIKE_NAME)].dropna()
	df_rgft = df_rgft[(df_rgft["prune_name"] != ROOTLIKE_NAME) & (df_rgft["rgft_name"] != ROOTLIKE_NAME)].dropna()
	df_prune.to_csv(outpath_prune)  # runover existing one to fill in all features
	df_rgft.to_csv(outpath_rgft)    # runover existing one to fill in all features
	
	return


def parse_neighbors_dirs(ds_path, outpath_prune, outpath_rgft, step_number, cp_internal=False, tree_type="bionj"):
	'''
	this function is used only when re-running SPR. otherwise, I need to parse the 'newick.csv' (and create truncated msas for subtrees)
	'''
	print("**** ", ds_path)

	msa_file = ds_path + MSA_PHYLIP_FILENAME
	all_trees = ds_path + REARRANGEMENTS_NAME + "s/"
	outpath_trees = TREES_PER_DS.format(ds_path, step_number)
	suf = "bionj" if tree_type == 'bionj' else 'br'  # if tree_type="random"


	res_dict_orig_tree = parse_phyml_stats_output(msa_file, ds_path + PHYML_STATS_FILENAME.format(suf))
	ll_orig_tree = float(res_dict_orig_tree["ll"])
	'''
	with open(ds_path + PHYML_STATS_FILENAME.format(suf)) as fpr:
		content = fpr.read()
		ll_orig_tree = float(re.search("Log-likelihood:\s+(.*)", content).group(1).strip())
	'''

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
					rearr_tree_path = SEP.join([tree_dirpath, "{}_phyml_{}_{}.txt".format(MSA_PHYLIP_FILENAME, "tree", OPT_TYPE)])
					cp_internal_names(rearr_tree_path, treepath_with_internal)


				ll_rearr = return_ll(tree_dirpath, OPT_TYPE)
				'''
				f = SEP.join([tree_dirpath, "{}_phyml_{}_{}.txt".format(MSA_PHYLIP_FILENAME, "stats", OPT_TYPE)])
				try:
					with open(f) as fpr:
						content = fpr.read()
						ll_rearr = float(re.search("Log-likelihood:\s+(.*)", content).group(1).strip())
				except:
					print(f)
					ll_rearr = None
				'''

				df.ix[ind, "prune_name"], df.ix[ind, "rgft_name"] = prune_name, rgft_name
				df.ix[ind, "orig_ds_ll"], df.ix[ind, "ll"] = ll_orig_tree, ll_rearr

	df.to_csv(outpath_prune)
	df.to_csv(outpath_rgft)
	df2.to_csv(outpath_trees)

	# to reduce inodes numbe r- delete subdirs after copying important content to 2 csvs:
	shutil.rmtree(all_trees, ignore_errors=True)

	return



def create_parsing_job(dataset_path, cp_internal, step_number, tree_type):
	job_name = "parse_features_ds"
	cmd = "python " + CODE_PATH + "summary/collect_SPR_features.py -ds " + dataset_path + " -st " + str(step_number) + " -ttype " + tree_type
	if cp_internal:
		cmd += " -cp "
	create_job_file.main(command=cmd, dirpath=dataset_path, sh_file=job_name + ".sh", multiply_jobs=False,
						 priority=-1, job_name=job_name)





if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='perform all SPR moves')
	parser.add_argument('--dataset_path', '-ds', default=None)
	parser.add_argument('--cp_internal', '-cp', default=False, action='store_true')
	parser.add_argument('--step_number', '-st', required=True)  # counting from 1
	parser.add_argument('--tree_type', '-ttype', default='bionj')  # could be bionj or random
	parser.add_argument('--index_to_start_run', '-istart', default=False)
	parser.add_argument('--nline_to_run', '-nlines', default=False)
	args = parser.parse_args()

	dataset_path = args.dataset_path
	if dataset_path:
		outpath_prune = SUMMARY_PER_DS.format(dataset_path, "prune", OPT_TYPE, args.step_number)
		outpath_rgft = SUMMARY_PER_DS.format(dataset_path, "rgft", "br", args.step_number)
	
		res = parse_neighbors_dirs(dataset_path, outpath_prune, outpath_rgft, args.step_number, args.cp_internal)
		collect_features(dataset_path, args.step_number, outpath_prune, outpath_rgft, args.tree_type)
	else:
		csv_path = SUMMARY_FILES_DIR + CHOSEN_DATASETS_FILENAME
		traverse_data_dirs(create_parsing_job, csv_path, (args.index_to_start_run, args.nline_to_run), args.cp_internal, args.step_number, args.tree_type)