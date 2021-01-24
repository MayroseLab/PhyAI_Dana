import sys

sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")


from defs import *
from my_utils.create_job_file import get_job_qsub_command
from summary.collect_SPR_features import *
from execute_programs.SPR_move import call_raxml_mem
from subprocess import Popen, PIPE, STDOUT
from parsing.parse_raxml_NG import parse_raxmlNG_content


EXAMPLE_DIRNAME = 'example3782/'



def run_raxml_mem_partitioned(tree_str, msa_tmpfile, log_filepath):
	RAXML_NG_SCRIPT = 'raxml-ng'

	# create tree file in memory and not in the storage:
	tree_rampath = "/dev/shm/" + str(random.random()) + str(random.random()) + "tree"  # the var is the str: tmp{dir_suffix}

	try:
		with open(tree_rampath, "w") as fpw:
			fpw.write(tree_str)

		p = Popen([RAXML_NG_SCRIPT, '--evaluate','--opt-branches', 'on',
				   '--opt-model', 'off', '--msa', msa_tmpfile, '--threads', '8', '--model', log_filepath, '--nofiles', '--redo', '--tree', tree_rampath],
				  stdout=PIPE, stdin=PIPE, stderr=STDOUT)
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


def index_ll_and_features(ds_path, outpath_prune, outpath_rgft, istart, nlines, NROWS):
	istart, nlines = int(istart)+1, int(nlines)
	skp_lst = [i for i in range(1, istart)] if not istart == 0 else []
	skp_lst2 = [i for i in range(istart+nlines, int(NROWS)+10)]
	skp_lst.extend(skp_lst2)
	dfr = pd.read_csv("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/data/training_datasets/{}newicks_step1.csv".format(EXAMPLE_DIRNAME),engine='python', index_col=0, skiprows=skp_lst)

	orig_ds_msa_file = ds_path + MSA_PHYLIP_FILENAME
	if not '59' in EXAMPLE_DIRNAME:
		stats_filepath = ds_path + PHYML_STATS_FILENAME.format('bionj')
		tree_file = ds_path + PHYML_TREE_FILENAME.format('bionj')
		params_dict = parse_phyml_stats_output(None, stats_filepath)
		freq, rates, pinv, alpha = [params_dict["fA"], params_dict["fC"], params_dict["fG"], params_dict["fT"]], [params_dict["subAC"], params_dict["subAG"], params_dict["subAT"], params_dict["subCG"],params_dict["subCT"], params_dict["subGT"]], params_dict["pInv"], params_dict["gamma"]
		orig_ds_ll = float(params_dict["ll"])
	else:
		tree_file = "/groups/itay_mayrose/danaazouri/PhyAI/DBset2/data/training_datasets/example59/tree_partitione_model_opt.txt"
		log_filepath = "/groups/itay_mayrose/danaazouri/PhyAI/DBset2/data/training_datasets/example59/partitioned/partitioned_output.raxml.bestModel"
		orig_ds_ll = -52052.674807


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

			if not '59' in EXAMPLE_DIRNAME:
				ll_rearr, rtime = call_raxml_mem(tree, orig_ds_msa_file, rates, pinv, alpha, freq)
			else:
				ll_rearr, rtime = run_raxml_mem_partitioned(tree, orig_ds_msa_file, log_filepath)

			df_prune.loc[ind, "time"], df_rgft.loc[ind, "time"] = rtime, rtime
			df_prune.loc[ind, "ll"], df_rgft.loc[ind, "ll"] = float(ll_rearr), float(ll_rearr)
			df_prune.loc[ind,"orig_ds_ll"], df_rgft.loc[ind,"orig_ds_ll"] = orig_ds_ll, orig_ds_ll

			features_restree_dict = calc_leaves_features(tree, "res", rgft_node_name=row["rgft_name"])
			df_prune = index_shared_features(df_prune, ind, row["prune_name"], "prune", features_prune_dicts_dict)
			df_rgft = index_shared_features(df_rgft, ind, row["rgft_name"], "rgft", features_rgft_dicts_dict)
			df_rgft = index_additional_rgft_features(df_rgft, ind, row["prune_name"], row["rgft_name"], features_restree_dict, features_prune_dicts_dict)  # also prune dict because for 2 features i didn't want to comp dict within each rgft iteration (needed to compute on the starting tree)

			df_rgft.loc[ind, FEATURES["res_bl"]] = features_restree_dict['res_bl']
			#df_rgft.loc[ind, FEATURES["res_tbl"]] = features_restree_dict['res_tbl']
			df_prune.to_csv(outpath_prune, mode='a', header=False)
			df_rgft.to_csv(outpath_rgft, mode='a', header=False)


	#df_prune = df_prune[(df_prune["prune_name"] != ROOTLIKE_NAME) & (df_prune["rgft_name"] != ROOTLIKE_NAME)]  # .dropna()
	#df_rgft = df_rgft[(df_rgft["prune_name"] != ROOTLIKE_NAME) & (df_rgft["rgft_name"] != ROOTLIKE_NAME)]  # .dropna()
	#df_prune.to_csv(outpath_prune)
	#df_rgft.to_csv(outpath_rgft)

	return


def submit_job_ll(istart, nlines, NROWS):
	print("**************************************   ", str(istart), str(nlines))
	job_name = "ll_large_dataset.sh"
	cmd = "python " + CODE_PATH + "my_utils/trying.py -istart " + str(istart) + " -nlines " + str(nlines) + " -nrows_total " + str(NROWS)

	qsub_cmd = get_job_qsub_command(job_name=job_name,
									command=cmd,
									error_files_path=DATA_PATH + EXAMPLE_DIRNAME +"error_files/",
									additional_params="-l pmem=15gb")
	os.system(qsub_cmd)



if __name__ == '__main__':
	df_train = pd.read_csv(SUMMARY_FILES_DIR + "play.csv", nrows=100)
	df_val = pd.read_csv(SUMMARY_FILES_DIR + "model_testing_validation_set.csv", nrows=100)
	#df_train = pd.read_csv("C:\\Users\\ItayMNB3\\Desktop\\play.csv", nrows=100)
	#df_val = pd.read_csv("C:\\Users\\ItayMNB3\\Desktop\\model_testing_validation_set.csv", nrows=100)
	from sklearn.ensemble import RandomForestRegressor
	
	def split_features_label(df, features):
		attributes_df = df[features].reset_index(drop=True)
		label_df = df['d_ll_merged'].reset_index(drop=True)
		
		x = np.array(attributes_df)
		y = np.array(label_df).ravel()
		
		return x, y
	
	features = FEATURES_MERGED
	features.remove(FEATURES['group_id'])
	X_train, y_train = split_features_label(df_train, features)
	X_test, y_test = split_features_label(df_val, features)
	
	model = RandomForestRegressor(oob_score=False, n_jobs=-1).fit(X_train,y_train)
	y_pred = model.predict(X_test)
	
	from statistics import mean
	def avg_sp_corr(y_true, y_pred):
		df_val["pred"] = np.array(y_pred)
		
		corrs = []
		grouped_df_by_ds = df_val.groupby(FEATURES["group_id"], sort=False)
		for group_id, df_by_ds in grouped_df_by_ds:
			temp_df = df_by_ds[['d_ll_merged', "pred"]]
			sp_corr = temp_df.corr(method='spearman').iloc[1, 0]
			corrs.append(sp_corr)
	
		return mean(corrs)
	
	mean_sp_corr = avg_sp_corr(None, y_pred)
	print(mean_sp_corr)
	
	###############
	from sklearn.model_selection import RandomizedSearchCV
	from sklearn.metrics import make_scorer
	sp_score = make_scorer(avg_sp_corr, greater_is_better=False)
	parameters = {'n_estimators': [25, 50, 70, 100],
				  'max_depth': [10, 30, 50, 70, None],
				  'min_samples_split': [2, 5, 10],
				  'min_samples_leaf': [1, 2, 5, 10],
				  'max_features': ["sqrt", "log2", 0.33, None],
				  'min_impurity_decrease': [0.0, 0.001, 0.01, 0.02, 0.1],
				  'ccp_alpha': [0.0, 0.1, 0.3, 0.5, 0.7]
				  }
	
	#grid_search = RandomizedSearchCV(model, parameters, cv=2, scoring=sp_score, n_jobs=-1)
	#opt_hp_model = grid_search.fit(X_train,y_train)
	#print(opt_hp_model.best_params_)
	
	#y_pred = opt_hp_model.predict(X_test)
	#mean_sp_corr = avg_sp_corr(None, y_pred)
	#print(mean_sp_corr)
	model2 = RandomForestRegressor(n_estimators=100, max_depth=70,
								  min_samples_split=10, min_samples_leaf=5,
								  max_features='sqrt',
								  min_impurity_decrease=0.001, ccp_alpha=0.1,
								  oob_score=False, n_jobs=-1).fit(X_train, y_train)
	y_pred = model2.predict(X_test)
	mean_sp_corr = avg_sp_corr(None, y_pred)
	print(mean_sp_corr)
	
	exit()
	
	
	parser = argparse.ArgumentParser(description='perform all SPR moves')
	parser.add_argument('--index_to_start_run', '-istart', default=False)
	parser.add_argument('--nline_to_run', '-nlines', default=False)
	parser.add_argument('--nrows_total_in_csv', '-nrows_total', default=False)
	args = parser.parse_args()

	if not args.index_to_start_run:
		df = pd.read_csv("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/data/training_datasets/{}newicks_step1.csv".format(EXAMPLE_DIRNAME)).dropna() #,index_col=0)

		NROWS = len(df)
		print(NROWS)
		print(df.tail())
		group_ids_full = df["group_id"]
		group_ids = group_ids_full.unique()
		for group in group_ids[1:]:
			s = df.index[df["group_id"] == group].tolist()
			submit_job_ll(s[0], len(s), NROWS)

	else:
		dataset_path = DATA_PATH + EXAMPLE_DIRNAME
		outpath_prune = SUMMARY_PER_DS.format(dataset_path + 'results_by_susbsets/', "prune", 'br', '1_subs_{}_{}'.format(args.index_to_start_run, args.nline_to_run))
		outpath_rgft = SUMMARY_PER_DS.format(dataset_path + 'results_by_susbsets/', "rgft", 'br', '1_subs_{}_{}'.format(args.index_to_start_run, args.nline_to_run))

		print(args.index_to_start_run, args.nline_to_run)
		index_ll_and_features(dataset_path, outpath_prune, outpath_rgft, args.index_to_start_run, args.nline_to_run, args.nrows_total_in_csv)


	'''
	dataset_path = DATA_PATH + EXAMPLE_DIRNAME
	df = pd.read_csv(dataset_path + "newicks_step1_with_ids.csv")  # -, index_col=1)

	group_ids_full = df["group_id"]
	group_ids = group_ids_full.unique()
	df_merged_prune, df_merged_rgft = pd.DataFrame(), pd.DataFrame()
	for group in group_ids:
		s = df.index[df["group_id"] == group].tolist()

		outpath_prune = pd.read_csv(SUMMARY_PER_DS.format(dataset_path + 'results_by_susbsets/', "prune", 'br','1_subs_{}_{}'.format(s[0], len(s))))
		outpath_rgft = pd.read_csv(SUMMARY_PER_DS.format(dataset_path + 'results_by_susbsets/', "rgft", 'br','1_subs_{}_{}'.format(s[0], len(s))))
		df_merged_prune = pd.concat([df_merged_prune, outpath_prune], ignore_index=True)
		df_merged_rgft = pd.concat([df_merged_rgft, outpath_rgft], ignore_index=True)
	df_merged_prune.to_csv(SUMMARY_PER_DS.format(dataset_path, "prune", 'br', '1'))
	df_merged_rgft.to_csv(SUMMARY_PER_DS.format(dataset_path, "rgft", 'br', '1'))
	exit()
	#'''
