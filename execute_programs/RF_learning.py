import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")
import warnings
warnings.filterwarnings("ignore")			# TEMP

from defs import *

from utils.tree_functions import get_total_branch_lengths
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
#from sklearn.preprocessing import *
from statistics import mean, median#from sklearn.metrics import *
#from sklearn.model_selection import train_test_split
#from sklearn.model_selection import cross_val_score
from figures.violin_for_grant import *
#from figures.confidence_interval_dts import plot_pred_true
#from figures.accXsize_boxplot import accXsize_boxplot
#from itertools import combinations
#import pickle
import joblib

pd.set_option('display.max_columns', 40)

ML_SOFTWARE_STATING_TREE = 'phyml'     # could be phyml | RAxML_NG
OPT_TYPE = "br"
KFOLD = 10     # "LOO"
GROUP_ID = 'group_id'
N_ESTIMATORS = 70
#MAX_DEPTH = 5
C = 95
FIGURES = False

FIRST_ON_RAND = False
FIRST_ON_SEC = False           # temp for running 1 on 2
FEATURE_SELECTION = False       # temp for running feature selection
SATURATION = False             # temp to asses saturation

N_DATASETS = 4200    # [1500,5858]


def score_rank(df_by_ds, sortby, locatein, random, scale_score):
	'''
	find the best tree in 'sortby' (e.g., predicted as the best) foreach dataset and locate its rank in 'locatein' (e.g., y_test)
	'''

	best_pred_ix = df_by_ds[sortby].idxmax()    # changed min to max!
	if random:
		best_pred_ix = np.random.choice(df_by_ds[sortby].index, 1, replace=False)[0]
	temp_df = df_by_ds.sort_values(by=locatein, ascending=False).reset_index()   # changed ascending to False
	best_pred_rank = min(temp_df.index[temp_df["index"] == best_pred_ix].tolist())
	best_pred_rank += 1  # convert from pythonic index to position
	
	if scale_score:
		best_pred_rank /= len(df_by_ds[sortby].index)   # scale the rank according to rankmax
		best_pred_rank *= 100


	return best_pred_rank


def save_df_ranks(df_by_ds, group_id, label):
	df_by_ds["ranked_label"] = df_by_ds[label].rank()
	df_by_ds["ranked_pred"] = df_by_ds["pred"].rank()

	cols = ["prune_name", "ranked_label", "ranked_pred", "sp_corr"]
	tmp = df_by_ds.sort_values(by="pred").reset_index()
	tmp["sp_corr"] = tmp[[label, "pred"]].corr(method='spearman').iloc[1,0]
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


def get_cumsun_preds(df_by_ds):
	df_by_ds["pred"] /= df_by_ds["pred"].sum()
	assert round(df_by_ds["pred"].sum()) == 1
	temp_df = df_by_ds.sort_values(by="pred", ascending=True).reset_index()      # ascending=True because the the sum() is negative so the order was changed
	sorted_preds = temp_df["pred"]
	cumsum_preds = sorted_preds.cumsum().values
	temp_df["pred"] = cumsum_preds
	
	return temp_df
	
	
def get_cumsum_threshold(df_by_ds, label):
	temp_df = get_cumsun_preds(df_by_ds)
	best_pred_ix = df_by_ds[label].idxmax()
	cumsum_true_best = temp_df[temp_df["index"] == best_pred_ix]["pred"].values[0]
	return cumsum_true_best


def calc_required_evaluations_score(grouped_df_by_ds, thresholds, c=C):
	cumsums = []
	threshold = np.percentile(thresholds, c)
	print("***",threshold)
	#plt.hist(thresholds)
	#plt.show()
	for group_id, df_by_ds in grouped_df_by_ds:
		cumulative_scores = get_cumsun_preds(df_by_ds)["pred"].values
		res = round(100 * (len(np.where(cumulative_scores < threshold)[0])) / len(cumulative_scores), 2)
		cumsums.append(res)

	return cumsums


def ds_scores(df, move_type, random, scale_score):
	rank_pred_by_ds, rank_test_by_ds = {}, {}

	label = LABEL.format(move_type)
	sp_corrs, r2s, errs_down, errs_up, all_true, all_preds, thresholds = [], [],[], [],[], [], []
	grouped_df_by_ds = df.groupby(FEATURES[GROUP_ID], sort=False)
	for group_id, df_by_ds in grouped_df_by_ds:
		if FIGURES:
			save_df_ranks(df_by_ds, group_id, label)

		rank_pred_by_ds[group_id] = score_rank(df_by_ds, "pred", label, random, scale_score)
		rank_test_by_ds[group_id] = score_rank(df_by_ds, label, "pred", random, scale_score)

		all_true.append(df_by_ds[label].mean())
		pred = df_by_ds["pred"].mean()
		all_preds.append(pred)
		r2s.append(0) #r2s.append(r2_score(df_by_ds[label], df_by_ds["pred"], multioutput='variance_weighted'))
		temp_df = df_by_ds[[label, "pred"]]
		sp_corr = temp_df.corr(method='spearman').iloc[1,0]
		if sp_corr:
			if sp_corr < 0:
				print(group_id, sp_corr)
			#sp_corrs.append(np.square(sp_corr))
			sp_corrs.append(sp_corr)
		else:
			sp_corrs.append(None)

		'''    # I don't use this score anymore
		# compute 'confidence score'
		cumsum_true_best = get_cumsum_threshold(df_by_ds, label)
		if cumsum_true_best:
			thresholds.append(cumsum_true_best)
		else:
			print("########", group_id)
	required_evaluations_scores = calc_required_evaluations_score(grouped_df_by_ds, thresholds, c=C)
	'''
	required_evaluations_scores = [0,0]    # I don't use this score anymore


	return rank_pred_by_ds, rank_test_by_ds, sp_corrs, r2s, required_evaluations_scores


def split_features_label(df, move_type, features):
	attributes_df = df[features].reset_index(drop=True)
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


def apply_RFR(df_test, df_train, move_type, features, cv=True):
	X_train, y_train = split_features_label(df_train, move_type, features)
	X_test, y_test = split_features_label(df_test, move_type, features)

	if not FEATURE_SELECTION and not SATURATION and not cv:
		model_path = SUMMARY_FILES_DIR + 'finalized_model_joblib.sav'
		if not os.path.exists(model_path):
			regressor = RandomForestRegressor(n_estimators=N_ESTIMATORS, max_features=0.33,  oob_score=True, n_jobs=-1).fit(X_train, y_train) # 0.33=nfeatures/3. this is like in R (instead of default=n_features)
			# save the model to disk
			joblib.dump(regressor, open(model_path, 'wb'))
		model = joblib.load(model_path)
	else:
		model = RandomForestRegressor(n_estimators=N_ESTIMATORS, max_features=0.33, oob_score=True, n_jobs=-1).fit(X_train, y_train)

	y_pred = model.predict(X_test)
	
	oob = model.oob_score_
	f_imp = model.feature_importances_

	all_DTs_pred = []


	return y_pred, all_DTs_pred, oob, f_imp


def truncate(df):
	groups_ids = df[FEATURES[GROUP_ID]].unique()
	#'''
	if DBSET == "2":
		n_other_dbs = 541
		selected_groups_ids = np.concatenate((np.random.choice(groups_ids[:-n_other_dbs], N_DATASETS-n_other_dbs, replace=False), groups_ids[-n_other_dbs:]))
		df = df[df[FEATURES[GROUP_ID]].isin(selected_groups_ids)]
		groups_ids = df[FEATURES[GROUP_ID]].unique()
	# '''
	kfold = len(groups_ids) if KFOLD=="LOO" else KFOLD
	assert len(groups_ids) >= kfold
	ndel = len(groups_ids) % kfold
	if ndel != 0:   # i removed datasets from the end, and not randomly. from some reason..
		for group_id in groups_ids[:-ndel-1:-1]:
			df = df[df[FEATURES[GROUP_ID]] != group_id]

	
	groups_ids = df[FEATURES[GROUP_ID]].unique()
	new_length = len(groups_ids) - ndel
	test_batch_size = int(new_length / kfold)

	return df.reset_index(drop=True), groups_ids, test_batch_size


def cross_validation_RF(df, move_type, features, trans=False, validation_set=None, random=False, scale_score=True):
	#'''
	df, groups_ids, test_batch_size = truncate(df)
	res_dict = {}
	oobs, f_imps, = [], []
	my_y_pred, imps = np.full(len(df), np.nan), np.full(len(df), np.nan)
	
	if not validation_set:

		for low_i in groups_ids[::test_batch_size]:
			low_i, = np.where(groups_ids == low_i)
			low_i = int(low_i)
			up_i = low_i + test_batch_size
	
			test_ixs = groups_ids[low_i:up_i]
			train_ixs = np.setdiff1d(groups_ids, test_ixs)
			df_test = df.loc[df[FEATURES[GROUP_ID]].isin(test_ixs)]
			df_train = df.loc[df[FEATURES[GROUP_ID]].isin(train_ixs)]
	
	
			y_pred, all_DTs_pred, oob, f_imp = apply_RFR(df_test, df_train, move_type, features, cv=True)
			
			oobs.append(oob)
			f_imps.append(f_imp)
			my_y_pred[df_test.index.values] = y_pred       # sort the predictions into a vector sorted according to the respective dataset
	
		df["pred"] = my_y_pred
	
	else:     # namely if validation set strategy, and not cross validation
		df_train = df
		if FIRST_ON_SEC:
			df_test = pd.read_csv(dirpath + LEARNING_DATA.format("all_moves", "2"))
		elif FIRST_ON_RAND:
			df_test = pd.read_csv(dirpath + LEARNING_DATA.format("all_moves", "1_random_starting"))
		else:   # a reg validation set
			df_test = pd.read_csv(dirpath + "model_testing_{}.csv".format(validation_set))

		df_test = fit_transformation(df_test, move_type, trans)  #.dropna()
		y_pred, all_DTs_pred, oob, f_imp = apply_RFR(df_test, df_train, move_type, features, cv=False)
		
		oobs.append(oob)
		f_imps.append(f_imp)
		df_test["pred"] = y_pred  # the predictions vec is the same lengths of test set
		df = df_test
	#'''
	
	
	rank_pred_by_ds, rank_test_by_ds, corrs, r2s, required_evaluations_scores = ds_scores(df, move_type, random, scale_score)

	# averaged over cv iterations
	res_dict['oob'] = sum(oobs) / len(oobs)
	res_dict['f_importance'] = sum(f_imps) / len(f_imps)
	# foreach dataset (namely arrays are of lengths len(sampled_datasets)
	res_dict["rank_first_pred"] = rank_pred_by_ds
	res_dict["rank_first_true"] = rank_test_by_ds
	res_dict["spearman_corr"] = corrs
	res_dict['%neighbors'] = required_evaluations_scores
	
	
	return res_dict, df


def fit_transformation(df, move_type, trans=False):
	# scale to the abs value of the starting tree ll
	df[LABEL.format(move_type)] /= -df["orig_ds_ll"]  # minus (abs) to preserve order. 'ascending' should be True in 'get_cumsun_preds' function

	if trans == 'rank':
		groups_ids = df[FEATURES[GROUP_ID]].unique()
		for group_id in groups_ids:
			df.loc[df[FEATURES[GROUP_ID]] == group_id, LABEL.format(move_type)] = df.loc[df[FEATURES[GROUP_ID]] == group_id, LABEL.format(move_type)].rank(ascending=True) #, pct=True)
	if trans == "standard":
		scaler = StandardScaler()
		scaler.fit(df[LABEL.format(move_type)].values.reshape(-1,1))
		df[LABEL.format(move_type)] = scaler.transform(df[LABEL.format(move_type)].values.reshape(-1,1)) + 100
	if trans == 'exp':
		df[LABEL.format(move_type)] = np.exp2(df[LABEL.format(move_type)]+1)
		#df[LABEL.format(move_type)] = df[LABEL.format(move_type)].transform(np.exp2)
	if trans == 'exp10':
		from scipy.special import exp10
		df[LABEL.format(move_type)] = exp10(df[LABEL.format(move_type)] + 1)

	return df


def parse_relevant_summaries_for_learning(df_orig, step_number, tree_type='bionj'):
	if 'example' in step_number:
		df_orig = pd.DataFrame(index=[0], columns=["path"])
		splt = step_number.split("_")
		df_orig.loc[0, "path"] = DATA_PATH + splt[0] + SEP
		step_number = "_".join(splt[1:]) if "_" in step_number else "1"

	ds_path_init = df_orig.iloc[0]["path"]
	cols1 = list(pd.read_csv(SUMMARY_PER_DS.format(ds_path_init, "prune", OPT_TYPE, step_number)))
	cols2 = list(pd.read_csv(SUMMARY_PER_DS.format(ds_path_init, "rgft", OPT_TYPE, step_number)))
	cols1.insert(1, "path")
	cols2.insert(1, "path")
	cols1.extend([FEATURES[GROUP_ID], FEATURES["group_tbl"]])  # add for group features
	cols2.extend([FEATURES[GROUP_ID], FEATURES["group_tbl"]])  # add for group features
	df1 = pd.DataFrame(index=np.arange(0), columns=cols1)
	df2 = pd.DataFrame(index=np.arange(0), columns=cols2)

	for i, row in df_orig.iterrows():
		ds_path = row["path"]
		ds_path = ds_path if tree_type == 'bionj' else ds_path + RANDOM_TREE_DIRNAME  # if == 'random
		suf = "bionj" if tree_type == 'bionj' else OPT_TYPE
		starting_tree = ds_path + PHYML_TREE_FILENAME.format(suf) if ML_SOFTWARE_STATING_TREE == 'phyml' else ds_path + RAXML_TREE_FILENAME
		if "pred" in step_number:
			starting_tree = ds_path + step_number + '.txt'
		if not os.path.exists(starting_tree):
			continue   # at least untill I finish debugging the failed ~150ds
		ds_tbl = get_total_branch_lengths(starting_tree)

		summary_per_ds1 = SUMMARY_PER_DS.format(ds_path, 'prune', OPT_TYPE, step_number)
		summary_per_ds2 = SUMMARY_PER_DS.format(ds_path, 'rgft', OPT_TYPE, step_number)
		print(summary_per_ds2)
		if os.path.exists(summary_per_ds2): # and FEATURES["bl"] in pd.read_csv(summary_per_ds2).columns:
			df_ds1 = pd.read_csv(summary_per_ds1)
			df_ds2 = pd.read_csv(summary_per_ds2)

			df_ds1.insert(1, "path", ds_path)
			df_ds2.insert(1, "path", ds_path)
			df_ds1[FEATURES[GROUP_ID]] = str(i)
			df_ds2[FEATURES[GROUP_ID]] = str(i)
			df_ds1[FEATURES["group_tbl"]] = ds_tbl
			df_ds2[FEATURES["group_tbl"]] = ds_tbl

			df1 = pd.concat([df1, df_ds1], ignore_index=True)
			df2 = pd.concat([df2, df_ds2], ignore_index=True)

	return df1, df2
	



def plot_cumulative_importance(f_list, importances, move_type):
	plt.figure()
	feature_importances = [(feature, round(importance, 4)) for feature, importance in zip(f_list, importances)]
	feature_importances = sorted(feature_importances, key=lambda x: x[1], reverse=True)
	x_values = list(range(len(importances.tolist())))
	sorted_importances = [importance[1] for importance in feature_importances]
	sorted_features = [importance[0] for importance in feature_importances]
	cumulative_importances = np.cumsum(sorted_importances)
	plt.plot(x_values, cumulative_importances, 'g-')
	plt.hlines(y=0.95, xmin=0, xmax=len(sorted_importances), color='r', linestyles='dashed')
	plt.xticks(x_values, sorted_features, rotation='vertical')
	plt.xlabel('Variable')
	plt.ylabel('Cumulative Importance')
	plt.title('Cumulative Importances')

	plt.tight_layout()
	plt.savefig(dirpath + "cumulative_importance_{}_{}.png".format(move_type, len(sorted_importances)))

	
	
def print_and_index_results(df_datasets, res_dict, move_type, features):
	
	#### score 1 ####
	spearman_corrs = res_dict['spearman_corr']
	df_datasets['corr'] = spearman_corrs
	print("\nsapearman corr:\n" + "mean:", mean([e for e in spearman_corrs if not math.isnan(e)]), ", median:",median(spearman_corrs))
	#print("\nsapearman corr:\n", spearman_corrs)
	
	#### score 2 + 3 ####
	res_vec1 = np.asarray(list(res_dict['rank_first_pred'].values())) if type(res_dict['rank_first_pred']) is dict else res_dict['rank_first_pred']
	res_vec2 = np.asarray(list(res_dict['rank_first_true'].values()))  if type(res_dict['rank_first_true']) is dict else res_dict['rank_first_true']
	res_vec1_scaled = res_vec1
	res_vec2_scaled = res_vec2
	df_datasets['best_predicted_ranking'] = res_vec1_scaled
	df_datasets['best_empirically_ranking'] = res_vec2_scaled
	print("\nbest predicted rank in true:\n" + "mean:",np.mean(res_vec1_scaled), ", median:", np.median(res_vec1_scaled))
	print("\nbest true rank in pred :\n" + "mean:",np.mean(res_vec2_scaled), ", median:", np.median(res_vec2_scaled))
	
	#### score 4 ####
	#required_evaluations = res_dict['%neighbors']
	#df_datasets['required_evaluations_0.95'] = required_evaluations
	#print("\nmean %neighbors (0.95): {}".format(sum(required_evaluations)/len(required_evaluations)))
	
	#'''### feature importance ####
	mean_importances = res_dict['f_importance']   # index in first row only (score foreach run and not foreach dataset)
	for i, f in enumerate(features):
		colname = "imp_" + f
		df_datasets.loc[0, colname] = mean_importances[i]
	#print("\nmean f importance:\n", np.column_stack((features, mean_importances)))
	plot_cumulative_importance(features, mean_importances, move_type)
	#'''
	#### additional information ####
	df_datasets.loc[0, 'oob'] = res_dict['oob']   # index in first row only (score foreach run and not foreach dataset)
	print("oob:", res_dict['oob'])
	print("ndatasets: ", len(res_vec1))
	
	print("##########################")
	return df_datasets


def sort_features(res_dict, features):
	feature_importances = [(feature, round(importance, 4)) for feature, importance in zip(features, res_dict['f_importance'])]
	feature_importances = sorted(feature_importances, key=lambda x: x[1], reverse=True)  # most important first
	sorted_importances = [importance[1] for importance in feature_importances]
	sorted_features = [importance[0] for importance in feature_importances]
	
	return sorted_importances, sorted_features
	
	
def extract_scores_dict(res_dict, df_with_scores):
	res_dict['rank_first_pred'], res_dict["rank_first_true"] = df_with_scores['best_predicted_ranking'].values, df_with_scores['best_empirically_ranking'].values
	res_dict['spearman_corr'], res_dict['oob'] = df_with_scores['corr'].values, df_with_scores.loc[0, 'oob']
	res_dict['f_importance'] = df_with_scores.loc[0, df_with_scores.columns[pd.Series(df_with_scores.columns).str.startswith('imp_')]].values

	return res_dict



if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='arrange data for learning and implement learning algo')
	parser.add_argument('--move_type', '-mt', default='prune')	 # could be 'prune' or 'rgft' or 'merged'
	parser.add_argument('--step_number', '-st', required=True) 	 # counting from 1   ## could also be 'example' !
	parser.add_argument('--validation_set', '-val', default=None) # could be None (default) | validation_set | example | example{n}
	parser.add_argument('--transform_target', '-trans', default=False)   # if trans, could be XX or YY
	parser.add_argument('--score_for_random', '-random', default=False, action='store_true')
	parser.add_argument('--scale_score', '-sscore', default=False, action='store_true')
	parser.add_argument('--tree_type', '-ttype', default='bionj')  # could be bionj or random
	args = parser.parse_args()

	dirpath = SUMMARY_FILES_DIR if platform.system() == 'Linux' else DATA_PATH
	df_orig = pd.read_csv(dirpath + CHOSEN_DATASETS_FILENAME, dtype=types_dict) #, nrows=30 OR skiprows=[i for i in range(200,5502)])

	move_type = args.move_type
	st = str(args.step_number)
	ifall = "all_moves"
	ifrandomstart = "" if args.tree_type == 'bionj' else "_random_starting"  # if == 'random'
	
	if not move_type == "merged":
		pass # no need anymore, antique analysis. if I need check versions before 27/9/20
		df_path = ''
	else:  # parse ALL neighbors to create a merged df off all features of all neighbors
		df_path = dirpath + LEARNING_DATA.format(ifall, st + ifrandomstart)

		if not os.path.exists(df_path):
			df_features_prune, df_features_rgft = parse_relevant_summaries_for_learning(df_orig, st, tree_type=args.tree_type)
			shared_cols = FEATURES_SHARED + ["path","prune_name","rgft_name","orig_ds_ll", "ll"]
			complete_df = df_features_prune.merge(df_features_rgft,on=shared_cols, left_index=True, right_index=True, suffixes=('_prune', '_rgft'))
			complete_df = complete_df.rename(columns={FEATURES[f]: FEATURES[f] + "_rgft" for f in FEATURES_RGFT_ONLY})
			complete_df[LABEL.format(move_type)] = complete_df[LABEL.format("prune")]
			complete_df.to_csv(df_path)
	if 'example' in st:   # if I only need to merge the summary file into model_testing_examplexx.csv
		exit()

	df_learning = pd.read_csv(df_path, dtype=types_dict)
	df_learning = fit_transformation(df_learning, move_type, trans=args.transform_target)

	features = FEATURES_PRUNE if move_type == "prune" else FEATURES_RGFT if move_type == "rgft" else FEATURES_MERGED
	features.remove(FEATURES[GROUP_ID])
	features_to_drop = []

	########################
	#'''
	val = args.validation_set
	suf = "_{st}_{valtype}".format(st=st, valtype=val) if val and not FIRST_ON_SEC else "_1st_on_2nd" if val else "_{}".format(st)
	ifsaturation = "" if not SATURATION else "_" + str(N_DATASETS)
	ifrank = "" if not args.transform_target else "_ytransformed_{}".format(args.transform_target)
	suf += ifsaturation + ifrank + ifrandomstart

	for i in range(len(features)):
		csv_with_scores = dirpath + SCORES_PER_DS.format(str(len(features))+ suf)
		csv_with_preds = dirpath + DATA_WITH_PREDS.format(str(len(features)) + suf)
		if not os.path.exists(csv_with_scores) or val:
			print("*@*@*@* scores for step{} with {} features are not available, thus applying learning".format(suf, len(features)))
			res_dict, df_out = cross_validation_RF(df_learning, move_type, features, trans=args.transform_target ,validation_set=val,random=args.score_for_random,scale_score=args.scale_score)
			df_out.to_csv(csv_with_preds)

			df_datasets = df_orig if not val else pd.read_csv(DIRPATH + "/validation_set2/summary_files/" + CHOSEN_DATASETS_FILENAME)
			df_datasets = df_datasets[df_datasets["path"].isin(df_out["path"].unique())]
		else:
			df_datasets = pd.read_csv(csv_with_scores)
			res_dict = extract_scores_dict({}, df_datasets)
		df_datasets = print_and_index_results(df_datasets, res_dict, move_type, features)
		df_datasets.to_csv(csv_with_scores)
		
		if not FEATURE_SELECTION:
			exit()
		else:
			sorted_importances, sorted_features = sort_features(res_dict, features)
			features = sorted_features[:-1]
			features_to_drop.append(sorted_features[-1])
			print("**dropped features:", features_to_drop)
	
	'''
	### temp - one-feature analysis
	for i,f in enumerate(features[2:]):
		print("***", f)
		csv_with_scores = dirpath + SCORES_PER_DS.format(f + '_' +str(i+1))
		if not os.path.exists(csv_with_scores):
			res_dict, df_out = cross_validation_RF(df_learning, move_type, [f], trans=args.transform_target, validation_set=None, random=args.score_for_random, scale_score=args.scale_score)
			df_out.to_csv(dirpath + DATA_WITH_PREDS.format(str(len(f)) + "_" + f))
			df_datasets = df_orig[df_orig["path"].isin(df_out["path"].unique())]
		else:
			df_datasets = pd.read_csv(csv_with_scores)
			res_dict = extract_scores_dict({}, df_datasets)
		df_datasets = print_and_index_results(df_datasets, res_dict, move_type, args.scale_score, [f])
		df_datasets.to_csv(csv_with_scores)
	#'''