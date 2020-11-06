import sys



sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")

from defs import *
from utils.create_job_file import get_job_qsub_command
from summary.collect_SPR_features import *
from execute_programs.SPR_move import call_raxml_mem



NROWS = 491100   #365400   #365380  ##491087


def index_ll_and_features(ds_path, outpath_prune, outpath_rgft, istart, nlines):
	istart, nlines = int(istart)+1, int(nlines)
	skp_lst = [i for i in range(1, istart)] if not istart == 0 else []
	skp_lst2 = [i for i in range(istart+nlines, NROWS)]
	skp_lst.extend(skp_lst2)
	dfr = pd.read_csv("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/data/training_datasets/example404/newicks_step1.csv", index_col=0, skiprows=skp_lst)

	orig_ds_msa_file = ds_path + MSA_PHYLIP_FILENAME
	stats_filepath = ds_path + PHYML_STATS_FILENAME.format('bionj')
	tree_file = ds_path + PHYML_TREE_FILENAME.format('bionj')
	params_dict = parse_phyml_stats_output(None, stats_filepath)
	freq, rates, pinv, alpha = [params_dict["fA"], params_dict["fC"], params_dict["fG"], params_dict["fT"]], [params_dict["subAC"], params_dict["subAG"], params_dict["subAT"], params_dict["subCG"],params_dict["subCT"], params_dict["subGT"]], params_dict["pInv"], params_dict["gamma"]
	orig_ds_ll = float(params_dict["ll"])

	features_prune_dicts_dict = calc_leaves_features(tree_file, "prune")
	df_prune, df_rgft = pd.DataFrame(), pd.DataFrame()
	for i, row in dfr.iterrows():
		ind = row.name
		tree = row["newick"]
		if row["rgft_name"] == SUBTREE2:  # namely the remaining subtree
			features_rgft_dicts_dict = calc_leaves_features(tree, "rgft")
		if not "subtree" in row["rgft_name"] and not ROOTLIKE_NAME in row["rgft_name"] and not ROOTLIKE_NAME in row["prune_name"]:
			prune_name, rgft_name = row["prune_name"], row['rgft_name']
			df_prune.loc[ind, "prune_name"], df_prune.loc[ind, "rgft_name"], df_rgft.loc[ind, "prune_name"], df_rgft.loc[ind, "rgft_name"] = prune_name, rgft_name, prune_name, rgft_name

			ll_rearr, rtime = call_raxml_mem(tree, orig_ds_msa_file, rates, pinv, alpha, freq)
			df_prune.loc[ind, "time"], df_rgft.loc[ind, "time"] = rtime, rtime
			df_prune.loc[ind, "ll"], df_rgft.loc[ind, "ll"] = float(ll_rearr), float(ll_rearr)
			df_prune.loc[ind,"orig_ds_ll"], df_rgft.loc[ind,"orig_ds_ll"] = orig_ds_ll, orig_ds_ll

			features_restree_dict = calc_leaves_features(tree, "res", rgft_node_name=row["rgft_name"])
			df_prune = index_shared_features(df_prune, ind, row["prune_name"], "prune", features_prune_dicts_dict)
			df_rgft = index_shared_features(df_rgft, ind, row["rgft_name"], "rgft", features_rgft_dicts_dict)
			df_rgft = index_additional_rgft_features(df_rgft, ind, row["prune_name"], row["rgft_name"], features_restree_dict, features_prune_dicts_dict)  # also prune dict because for 2 features i didn't want to comp dict within each rgft iteration (needed to compute on the starting tree)

			df_rgft.loc[ind, FEATURES["res_bl"]] = features_restree_dict['res_bl']
			df_rgft.loc[ind, FEATURES["res_tbl"]] = features_restree_dict['res_tbl']

	df_prune = df_prune[(df_prune["prune_name"] != ROOTLIKE_NAME) & (df_prune["rgft_name"] != ROOTLIKE_NAME)]  # .dropna()
	df_rgft = df_rgft[(df_rgft["prune_name"] != ROOTLIKE_NAME) & (df_rgft["rgft_name"] != ROOTLIKE_NAME)]  # .dropna()
	df_prune.to_csv(outpath_prune)  # runover existing one to fill in all features
	df_rgft.to_csv(outpath_rgft)  # runover existing one to fill in all features

	return


def submit_job_ll(istart, nlines):
	print("**************************************   ", str(istart), str(nlines))
	job_name = "index_ll_large_dataset.sh"
	cmd = "python " + CODE_PATH + "utils/trying.py -istart " + str(istart) + " -nlines " + str(nlines)

	qsub_cmd = get_job_qsub_command(job_name=job_name,
									command=cmd,
									error_files_path=DATA_PATH + "example404/error_files/")
	os.system(qsub_cmd)



if __name__ == '__main__':
	#df = pd.read_csv("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/summary_files/with_preds_merged_20_1_ml_minus1_set_4200_ytransformed_exp.csv")
	#temp_df = df.sort_values(by='pred', ascending=False).reset_index()
	#print(temp_df['d_ll_prune'].head(30), "\n")
	#temp_df = df.sort_values(by='d_ll_merged', ascending=False).reset_index()
	#print(temp_df['d_ll_prune'].head(30), "\n")

	parser = argparse.ArgumentParser(description='perform all SPR moves')
	parser.add_argument('--index_to_start_run', '-istart', default=False)
	parser.add_argument('--nline_to_run', '-nlines', default=False)
	args = parser.parse_args()

	if not args.index_to_start_run:
		df = pd.read_csv("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/data/training_datasets/example404/newicks_step1_with_ids.csv")

		group_ids_full = df["group_id"]
		group_ids = group_ids_full.unique()
		for group in group_ids[9:]:
			s = df.index[df["group_id"] == group].tolist()
			submit_job_ll(s[0], len(s))
	else:
		dataset_path = DATA_PATH + 'example404/'
		outpath_prune = SUMMARY_PER_DS.format(dataset_path + 'results_by_susbsets/', "prune", 'br', '1_subs_{}_{}'.format(args.index_to_start_run, args.nline_to_run))
		outpath_rgft = SUMMARY_PER_DS.format(dataset_path + 'results_by_susbsets/', "rgft", 'br', '1_subs_{}_{}'.format(args.index_to_start_run, args.nline_to_run))

		print(args.index_to_start_run, args.nline_to_run)
		index_ll_and_features(dataset_path, outpath_prune, outpath_rgft, args.index_to_start_run, args.nline_to_run)
