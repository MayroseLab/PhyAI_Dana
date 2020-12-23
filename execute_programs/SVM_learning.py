import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")
import warnings
warnings.filterwarnings("ignore")			# TEMP


from defs import *

from my_utils.tree_functions import get_total_branch_lengths
from sklearn import svm
from sklearn.base import RegressorMixin
from sklearn.metrics import r2_score
from execute_programs.RF_learning import *
from sklearn import linear_model
from sklearn.neighbors import KNeighborsRegressor
from sklearn.neural_network import MLPRegressor


GROUP_ID = 'group_id'
SATURATION = False



def apply_SVR(df_test, df_train, move_type, features):
	X_train, y_train = split_features_label(df_train, move_type, features)
	X_test, y_test = split_features_label(df_test, move_type, features)

	#clf = svm.SVR(epsilon=0.000005)    # worse, but not so bad, takes ages
	clf = linear_model.BayesianRidge()	# almost comparable to RandomForest! (maybe better on smaller trees)   **READ**
	#clf = linear_model.SGDRegressor()	# worse (fast)
	#clf = linear_model.Lasso()			# worse (fast)
	#clf = KNeighborsRegressor()		# worse
	###clf = MLPRegressor()	# WRONG		# bad as random (NN Multi-layer Perceptron regressor) 	# WRONG
	clf.fit(X_train, y_train)
	y_pred = clf.predict(X_test)

	return y_pred



def cross_validation_SVM(df, move_type, features, random, scale_score):
	df, groups_ids, test_batch_size = truncate(df)
	
	rank_score_by_ds, res_dict = {}, {}
	my_y_pred = np.full(len(df), np.nan)
	for low_i in groups_ids[::test_batch_size]:
		low_i, = np.where(groups_ids == low_i)
		low_i = int(low_i)
		up_i = low_i + test_batch_size

		test_ixs = groups_ids[low_i:up_i]
		train_ixs = np.setdiff1d(groups_ids, test_ixs)
		df_test = df.loc[df[FEATURES[GROUP_ID]].isin(test_ixs)]
		df_train = df.loc[df[FEATURES[GROUP_ID]].isin(train_ixs)]
		y_pred = apply_SVR(df_test, df_train, move_type, features)

		my_y_pred[df_test.index.values] = y_pred       # sort the predictions into a vector sorted according to the respective dataset
		
	df["pred"] = my_y_pred
	rank_pred_by_ds, rank_test_by_ds, corrs, r2s, required_evaluations_scores = ds_scores(df, move_type, random, scale_score)


	# foreach dataset (namely arrays are of lengths len(sampled_datasets)
	res_dict["rank_first_pred"] = rank_pred_by_ds
	res_dict["rank_first_true"] = rank_test_by_ds
	res_dict["spearman_corr"] = corrs
	res_dict['%neighbors'] = required_evaluations_scores

	return res_dict, df



def print_and_index_results_SVM(df_datasets, res_dict):
	#### score 1 ####
	spearman_corrs = res_dict['spearman_corr']
	df_datasets['corr'] = spearman_corrs
	print("\nsapearman corr:\n" + "mean:", mean([e for e in spearman_corrs if not math.isnan(e)]), ", median:",
	      median(spearman_corrs))
	
	#### score 2 + 3 ####
	res_vec1 = np.asarray(list(res_dict['rank_first_pred'].values())) if type(res_dict['rank_first_pred']) is dict else \
	res_dict['rank_first_pred']
	res_vec2 = np.asarray(list(res_dict['rank_first_true'].values())) if type(res_dict['rank_first_true']) is dict else \
	res_dict['rank_first_true']
	scores_range = (1, 100)  # for MinMaxScaler
	res_vec1_scaled = ((res_vec1 - res_vec1.min(axis=0)) / (res_vec1.max(axis=0) - res_vec1.min(axis=0))) * (
				scores_range[1] - scores_range[0]) + scores_range[0]
	res_vec2_scaled = ((res_vec2 - res_vec2.min(axis=0)) / (res_vec2.max(axis=0) - res_vec2.min(axis=0))) * (
				scores_range[1] - scores_range[0]) + scores_range[0]
	df_datasets['best_predicted_ranking'] = res_vec1_scaled
	df_datasets['best_empirically_ranking'] = res_vec2_scaled
	print("\nbest predicted rank in true:\n" + "mean:", np.mean(res_vec1_scaled), ", median:",
	      np.median(res_vec1_scaled))
	print("\nbest true rank in pred :\n" + "mean:", np.mean(res_vec2_scaled), ", median:", np.median(res_vec2_scaled))
	
	#### score 4 ####
	res_vec2_scaled.sort()
	sorted_true_scaled = res_vec2_scaled[::-1]
	index95 = int(0.05 * len(sorted_true_scaled) - 1)  # index for loc 0.05 = 95%
	required_evaluations = res_dict['%neighbors']
	df_datasets['required_evaluations_0.95'] = required_evaluations
	print("\nIn 0.95: {}%".format(sorted_true_scaled[index95]))
	print("\nmean %neighbors (0.95): {}".format(sum(required_evaluations) / len(required_evaluations)))

	print("ndatasets: ", len(res_vec1))
	print("##########################")
	
	return df_datasets




if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='arrange data for learning and implement learning algo')
	parser.add_argument('--move_type', '-mt', default=False)	 # could be 'prune' or 'rgft'
	parser.add_argument('--step_number', '-st', required=True)  # counting from 1
	parser.add_argument('--score_for_random', '-random', default=False, action='store_true')
	parser.add_argument('--scale_score', '-sscore', default=False, action='store_true')
	parser.add_argument('--transform_target', '-trans', default=False)  # if trans, could be XX or YY
	parser.add_argument('--tree_type', '-ttype', default='bionj')  # could be bionj or random
	args = parser.parse_args()

	move_type = args.move_type
	dirpath = SUMMARY_FILES_DIR if platform.system() == 'Linux' else DATA_PATH
	df_orig = pd.read_csv(dirpath + CHOSEN_DATASETS_FILENAME, dtype=types_dict)
	df_orig = df_orig[df_orig["path"] != "/groups/itay_mayrose/danaazouri/PhyAI/DBset2/data/training_datasets/balibase_RV40_BB40049/"]
	
	
	suf = "_BayesianRidgeRegression_{}".format(args.step_number)
	ifsaturation = "" if not SATURATION else "_" + str(N_DATASETS)
	ifrank = "" if not args.transform_target else "_ytransformed_{}".format(args.transform_target)
	ifrandomstart = "" if args.tree_type == 'bionj' else "_random_starting"  # if == 'random'
	suf += ifsaturation + ifrank + ifrandomstart
	
	df_path = dirpath + LEARNING_DATA.format("all_moves", args.step_number + ifrandomstart)
	df_learning = pd.read_csv(df_path)
	df_learning = fit_transformation(df_learning, move_type, trans=False)
	
	features = FEATURES_PRUNE if move_type == "prune" else FEATURES_RGFT if move_type == "rgft" else FEATURES_MERGED
	features.remove(FEATURES[GROUP_ID])
	features_to_drop = []
	
	csv_with_scores = dirpath + SCORES_PER_DS.format(str(len(features)) + suf)
	if not os.path.exists(csv_with_scores):
		print("*@*@*@* scores for step{} with {} features are not available, thus applying learning".format(suf,len(features)))
		res_dict, df_out = cross_validation_SVM(df_learning, move_type, features, random=args.score_for_random, scale_score=args.scale_score)
		df_out.to_csv(dirpath + DATA_WITH_PREDS.format(str(len(features)) + suf))
		df_datasets = df_orig
		df_datasets = df_datasets[df_datasets["path"].isin(df_out["path"].unique())]
		
		df_datasets = print_and_index_results_SVM(df_datasets, res_dict)
		df_datasets.to_csv(csv_with_scores)
	else:
		print("*@*@*@* scores for step{} with {} features already exist".format(suf,len(features)))

