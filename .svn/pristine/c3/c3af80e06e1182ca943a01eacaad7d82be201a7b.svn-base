import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")
import warnings
warnings.filterwarnings("ignore")			# TEMP


from defs import *

from utils.tree_functions import get_total_branch_lengths
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn import svm
from sklearn.metrics import *
from sklearn import preprocessing
from sklearn.model_selection import cross_val_score
from statistics import mean, median
from figures.violin_for_grant import *
from figures.confidence_interval_dts import plot_pred_true
from figures.accXsize_boxplot import accXsize_boxplot


OPT_TYPE = "br"
KFOLD = 10     # "LOO"
GROUP_ID = 'group_id'
N_ESTIMATORS = 100
FIGURES = False


def score_rank(df_by_ds, sortby, locatein):
	'''
	find the best tree in 'sortby' (e.g., predicted as the best) foreach dataset and locate its rank in 'locatein' (e.g., y_test)
	'''
	# to do - generalize var names - not only the first direction
	best_pred_ix = df_by_ds[sortby].idxmin()
	temp_df = df_by_ds.sort_values(by=locatein).reset_index()
	best_pred_rank = min(temp_df.index[temp_df["index"] == best_pred_ix].tolist())
	best_pred_rank += 1  # convert from pythonic index to position

	return best_pred_rank


def save_df_ranks(df_by_ds, group_id, label):
	df_by_ds["ranked_label"] = df_by_ds[label].rank()
	df_by_ds["ranked_pred"] = df_by_ds["pred"].rank()

	cols = ["prune_name", "ranked_label", "ranked_pred", "sp_corr"]
	tmp = df_by_ds.sort_values(by="pred").reset_index()
	tmp["sp_corr"] = tmp[[label, "pred"]].corr(method='spearman').ix[1,0]
	tmp[cols].to_csv(DATA_PATH + str(group_id) + ".csv")


def dts_assessment(df_by_ds, label, pred, errs_down, errs_up, ranked=False):
	dts_scores = []
	for k in range(N_ESTIMATORS):
		col = "dt" + str(k)
		if not ranked:
			dts_scores.append(df_by_ds[col].mean())
		else:
			dts_scores.append(score_rank(df_by_ds, col, label))
			pred = 1
	res, err_down, err_up = confidence_score(dts_scores, pred)
	errs_down.append(err_down)
	errs_up.append(err_up)

	return errs_down, errs_up


def ds_scores(df, move_type):
	rank_pred_by_ds, rank_test_by_ds = {}, {}

	label = LABEL.format(move_type)
	sp_corrs, r2s, errs_down, errs_up, all_true, all_preds = [], [],[], [],[], []
	grouped_df_by_ds = df.groupby(FEATURES[GROUP_ID], sort=False)
	for group_id, df_by_ds in grouped_df_by_ds:
		if FIGURES:
			save_df_ranks(df_by_ds, group_id, label)

		rank_pred_by_ds[group_id] = score_rank(df_by_ds, "pred", label)
		rank_test_by_ds[group_id] = score_rank(df_by_ds, label, "pred")

		all_true.append(df_by_ds[label].mean())
		pred = df_by_ds["pred"].mean()
		all_preds.append(pred)


		r2s.append(r2_score(df_by_ds[label], df_by_ds["pred"], multioutput='variance_weighted'))
		temp_df = df_by_ds[[label, "pred"]]
		sp_corrs.append(temp_df.corr(method='spearman').ix[1,0])

	if FIGURES:
		ranked = False
		errs_down, errs_up = dts_assessment(df_by_ds, label, pred, errs_down, errs_up, ranked)
		X = list(rank_pred_by_ds.keys())
		if ranked:
			plot_pred_true(X, list(rank_pred_by_ds.values()), np.ones(len(X)), errs_down, errs_up)
		else:
			plot_pred_true(X, all_preds, all_true, errs_down, errs_up)


	return rank_pred_by_ds, rank_test_by_ds, sp_corrs, r2s


def split_features_label(df, move_type, feature_to_drop=GROUP_ID):
	features = FEATURES_PRUNE if move_type == "prune" else FEATURES_RGFT
	attributes_df = df[features].reset_index(drop=True)
	attributes_df = attributes_df.drop(FEATURES[feature_to_drop], axis=1)  # group_feature can be 'group_tbl' or 'group_id'
	label_df = df[LABEL.format(move_type)].reset_index(drop=True)

	x = np.array(attributes_df)
	y = np.array(label_df).ravel()

	return x, y


def confidence_score(all_samples, ypred, percentile=90):
	err_down = np.percentile(all_samples, (100 - percentile) / 2.0)
	err_up = np.percentile(all_samples,  100 - (100 - percentile) / 2.0)
	#std_vec = err_up-err_down/ypred
	#return std_vec#return std_vec
	if err_down <= ypred <= err_up:
		return 1, err_down, err_up
	return 0, err_down, err_up


def apply_RFR(df_test, df_train, move_type):
	X_train, y_train = split_features_label(df_train, move_type)
	X_test, y_test = split_features_label(df_test, move_type)

	regressor = RandomForestRegressor(n_estimators=N_ESTIMATORS, max_features=0.33,  oob_score=True).fit(X_train, y_train) # 0.33=nfeatures/3. this is like in R (instead of default=n_features)
	y_pred = regressor.predict(X_test)
	oob = regressor.oob_score_
	f_imp = regressor.feature_importances_

	all_DTs_pred = [t.predict(X_test) for t in regressor.estimators_]
	#dev_vec = confidence_score(all_DTs_pred, y_pred, percentile=90)


	return y_pred, all_DTs_pred, oob, f_imp


def truncate(df):
	groups_ids = df[FEATURES[GROUP_ID]].unique()
	kfold = len(groups_ids) if KFOLD=="LOO" else KFOLD
	assert len(groups_ids) >= kfold
	ndel = len(groups_ids) % kfold
	new_length = len(groups_ids) - ndel
	test_batch_size = int(new_length / kfold)
	if ndel != 0:
		for group_id in groups_ids[:-ndel-1:-1]:
			df = df[df[FEATURES[GROUP_ID]] != group_id]
		groups_ids = df[FEATURES[GROUP_ID]].unique()

	return df, groups_ids, test_batch_size


def cross_validation_RF(df, move_type):
	df, groups_ids, test_batch_size = truncate(df)

	oobs, f_imps, = [], []
	res_dict = {}
	my_y_pred = np.full(len(df), np.nan)
	all_dts = [np.full(len(df), np.nan) for i in range(N_ESTIMATORS)]

	k = 0.5
	for low_i in groups_ids[::test_batch_size]:
		low_i, = np.where(groups_ids == low_i)
		low_i = int(low_i)
		up_i = low_i + test_batch_size

		test_ixs = groups_ids[low_i:up_i]
		train_ixs = np.setdiff1d(groups_ids, test_ixs)
		df_test = df.loc[df[FEATURES[GROUP_ID]].isin(test_ixs)]
		df_train = df.loc[df[FEATURES[GROUP_ID]].isin(train_ixs)]

		#df_test = df_test.sample(frac=k)
		#df_train = df_train.sample(frac=k)
		#df_train = df_train[df_train["ntaxa"] == 60]

		y_pred, all_DTs_pred, oob, f_imp = apply_RFR(df_test, df_train, move_type)
		oobs.append(oob)
		f_imps.append(f_imp)
		my_y_pred[df_test.index.values] = y_pred       # sort the predictions into a vector sorted according to the respective dataset
		for ix in range(len((all_dts))):
			all_dts[ix][df_test.index.values] = all_DTs_pred[ix]


	for k, dt_preds in enumerate(all_dts):
		df["dt" + str(k)] = dt_preds
	df["pred"] = my_y_pred
	rank_pred_by_ds, rank_test_by_ds, corrs, r2s = ds_scores(df, move_type)

	res_dict["rank_first_pred"] = rank_pred_by_ds
	res_dict["rank_first_true"] = rank_test_by_ds
	res_dict["spearman_corr"] = corrs
	res_dict['mean oob'] = sum(oobs) / len(oobs)
	res_dict['mean f importance'] = sum(f_imps) / len(f_imps)

	return res_dict



def fit_transform(df, move_type, rank=False):
	groups_ids = df[FEATURES[GROUP_ID]].unique()
	#print(len(groups_ids))
	#print(groups_ids)
	#exit()
	for group_id in groups_ids:
		#scaling_factor = df[df[FEATURES[GROUP_ID]] == group_id][LABEL.format(move_type)].min()
		scaling_factor = df[df[FEATURES[GROUP_ID]] == group_id]["orig_ds_ll"].iloc[0]
		df.loc[df[FEATURES[GROUP_ID]] == group_id, LABEL.format(move_type)] /= scaling_factor
		if rank:
			df.loc[df[FEATURES[GROUP_ID]] == group_id, LABEL.format(move_type)] = df.loc[df[FEATURES[GROUP_ID]] == group_id, LABEL.format(move_type)].rank()

	end = int(len(df)*2/3)
	return df




def parse_relevant_summaries_for_learning(df_orig, outpath, move_type, step_number, all_moves=False):
	for i, row in df_orig.iterrows():
		ds_path = row["path"]
		ds_tbl = get_total_branch_lengths(ds_path + PHYML_TREE_FILENAME.format('bionj'))
		summary_per_ds = SUMMARY_PER_DS.format(ds_path, move_type, OPT_TYPE, step_number)
		print(summary_per_ds)
		if os.path.exists(summary_per_ds) and FEATURES["bl"] in pd.read_csv(summary_per_ds).columns:
			df_ds = pd.read_csv(summary_per_ds)
			if i == 0:
				cols = list(df_ds)
				if all_moves:
					cols = [col+"_"+move_type for col in cols]
				cols.insert(1, "path")
				cols.extend([FEATURES[GROUP_ID], FEATURES["group_tbl"]])		# add for group features
				df = pd.DataFrame(index=np.arange(0), columns=cols)

			grouped = df_ds.groupby("{}_name".format(move_type), sort=False)
			for j, (name, group) in enumerate(grouped):
				if all_moves:
					for k, (name_k, row_k) in enumerate(group.iterrows()):
						index = str(i) + "," + str(j+k)
						row_to_write = list(row_k.values)
						row_to_write.insert(1, ds_path)
						row_to_write.extend([str(i), ds_tbl])
						df.ix[index] = row_to_write
				else:
					best_row_group = list(group.ix[group[LABEL.format(move_type)].astype(float).idxmin()].values)
					best_row_group.insert(1, ds_path)
					best_row_group.extend([str(i), ds_tbl])							# add group features
					df.ix[str(i) + "," + str(j)] = best_row_group
		else:
			print(summary_per_ds, "does not exist!")
	df.to_csv(outpath)



def print_results(res_dict):
	print("\napearman corr:\n" + "mean:", mean(res_dict['spearman_corr']), ", median:", median(res_dict['spearman_corr']))
	print("\nmean oob:\n", res_dict['mean oob'])
	f_list = FEATURES_PRUNE if move_type == "prune" else FEATURES_RGFT
	f_list.remove(FEATURES[GROUP_ID])
	print("\nmean f importance:\n", np.column_stack((f_list, res_dict['mean f importance'])))

	res_vec1 = np.asarray(list(res_dict['rank_first_pred'].values()))
	res_vec2 = np.asarray(list(res_dict['rank_first_true'].values()))

	print("\nbest predicted rank in true:\n" + "mean:",np.mean(res_vec1), ", median:", np.median(res_vec1))
	print("\nbest true rank in pred :\n" + "mean:",np.mean(res_vec2), ", median:",  np.median(res_vec2))

	res_vec2.sort()
	sorted_true = res_vec2[::-1]
	index = int(0.01*len(sorted_true)-1)  # index for loc 0.01 = 99%
	print("\nPrecentages of samples needed: {}%".format(sorted_true[index]))
	# correct - this is a precentage only if i scale the rank with ntaxa

	if FIGURES:
		plot_violin({"rank of best predicted\nSPR move": res_vec1}, {"predicted rank of best\nSPR move": res_vec2})
		plot_hist({"rank of best predicted\nSPR move": res_vec1}, {"predicted rank of best\nSPR move": res_vec2})
		accXsize_boxplot(res_vec1)

	return


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='arrange data for learning and implement learning algo')
	parser.add_argument('--move_type', '-mt', default='prune')	 # could be 'prune' or 'rgft'
	parser.add_argument('--step_number', '-st', required=True)	 # counting from 1
	parser.add_argument('--all_moves', '-all', default=False, action='store_true')
	args = parser.parse_args()

	move_type = args.move_type
	dirpath = SUMMARY_FILES_DIR if platform.system() == 'Linux' else DATA_PATH

	df_orig = pd.read_csv(dirpath + CHOSEN_DATASETS_FILENAME)
	if not args.all_moves:
		df_path = dirpath + LEARNING_DATA.format(move_type, str(args.step_number))
		if not os.path.exists(df_path):
			parse_relevant_summaries_for_learning(df_orig, df_path, move_type, args.step_number)
	else:
		df_path = dirpath + LEARNING_DATA.format("all_moves", str(args.step_number))
		df_prune_features = dirpath + LEARNING_DATA.format("prune", str(args.step_number))
		df_rgft_features = dirpath + LEARNING_DATA.format("rgft", str(args.step_number))
		print("x")
		'''
		if not os.path.exists(df_path):
			# TODO: add move_type to features names.   DONE
			parse_relevant_summaries_for_learning(df_orig, df_prune_features, "prune", args.step_number, all_moves=True)
			parse_relevant_summaries_for_learning(df_orig, df_rgft_features, "rgft", args.step_number, all_moves=True)
			# TODO: df_path = pd.concat([df_prune_features, df_rgft_features], axis=1, sort=False).  OR:
			df_prune_features.merge(df_rgft_features,on=[, ], left_index=True, right_index=True, suffixes=('_prune', '_rgft'))
		'''

	df = pd.read_csv(df_path)
	#df = df[df["ntaxa"].astype('int64') != 60].reset_index(drop=True)

	df = fit_transform(df, move_type) #, rank=True)
	res_dict = cross_validation_RF(df, move_type)

	print_results(res_dict)