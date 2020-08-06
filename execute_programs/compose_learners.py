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
from execute_programs.RF_learning import score_rank

GROUP_ID = 'group_id'


def scale_to_prob(df, group_id):
	scaling_factor = df[df[FEATURES[GROUP_ID]] == group_id]["pred"].sum()
	df.loc[df[FEATURES[GROUP_ID]] == group_id, "pred"] /= scaling_factor
	#assert round(df.loc[df[FEATURES[GROUP_ID]] == group_id, "pred"].sum()) == 1
	if not round(df.loc[df[FEATURES[GROUP_ID]] == group_id, "pred"].sum()) == 1:
		print(group_id, df.loc[df[FEATURES[GROUP_ID]] == group_id, "pred"].sum())

	return df


def preds_to_probs(df_prune, df_rgft):
	'''
	convert ll predicted values to probs
	'''
	groups_ids = df_prune[FEATURES[GROUP_ID]].unique()
	for group_id in groups_ids:
		df_prune = scale_to_prob(df_prune, group_id)
		df_rgft = scale_to_prob(df_rgft, group_id)

	return df_prune, df_rgft


def compose_learners(df_prune, df_rgft, operator):
	#df_prune = df_prune.set_index('Unnamed: 0.1.1')
	#print(df_prune["pred"].astype(float).idxmax())
	#pd.set_option('display.max_columns', None)

	df_comp = pd.DataFrame(index=np.arange(0), columns=[])
	for i, row in df_prune.iterrows():
		ind_col = 'Unnamed: 0.1.1'
		ind = row[ind_col]
		path = row["path"]
		if len(df_rgft.loc[(df_rgft[ind_col] == ind) & (df_rgft["path"] == path) ,"pred"].values) > 0:
			for i,val in enumerate(df_rgft.loc[(df_rgft["path"] == path) & (df_rgft[ind_col] == ind), "pred"].values):
				df_comp.ix[(ind + "," +str(i)), "path"] = path
				df_comp.ix[(ind + "," +str(i)), FEATURES[GROUP_ID]] = row[FEATURES[GROUP_ID]]
				df_comp.ix[(ind + "," +str(i)), LABEL.format("prune")] = row[LABEL.format("prune")]

				if operator == "mult":
					df_comp.ix[(ind + "," + str(i)), "composed_pred"] = row["pred"] * (
					df_rgft.loc[(df_rgft["path"] == path) & (df_rgft[ind_col] == ind), "pred"].values[i])
				elif operator == "dev":
					df_comp.ix[(ind + "," + str(i)), "composed_pred"] = row["pred"] / (
					df_rgft.loc[(df_rgft["path"] == path) & (df_rgft[ind_col] == ind), "pred"].values[i])
				elif operator == "add":
					df_comp.ix[(ind + "," + str(i)), "composed_pred"] = row["pred"] + (
					df_rgft.loc[(df_rgft["path"] == path) & (df_rgft[ind_col] == ind), "pred"].values[i])
				else:
					print("no defined function for {} operator".format(operator))

	print(df_comp.head())
	df_comp.to_csv(dirpath + "comp_preds.csv")

	return df_comp.dropna()   # adding to 'df_prune' is just a random choice


def performance_evaluation(df_comp):
	label = LABEL.format("prune")
	rank_pred_by_ds, rank_test_by_ds = {}, {}
	grouped_df_by_ds = df_comp.groupby(FEATURES[GROUP_ID], sort=False)
	
	corrs = []
	for group_id, df_by_ds in grouped_df_by_ds:
		rank_pred_by_ds[group_id] = score_rank(df_by_ds, "composed_pred", label)
		rank_test_by_ds[group_id] = score_rank(df_by_ds, label, "composed_pred")
		
		corrs.append(df_by_ds[[label, "composed_pred"]].corr(method='spearman').ix[1,0])
	print(mean(corrs))
	print(median(corrs))
	return rank_pred_by_ds, rank_test_by_ds


def print_results(rank_pred_by_ds, rank_test_by_ds):
	print(len(rank_pred_by_ds))
	# MinMax scaler
	range = (1, 100)
	res_vec1 = np.asarray(list(rank_pred_by_ds.values()))
	res_vec2 = np.asarray(list(rank_test_by_ds.values()))
	res_vec1_scaled = ((res_vec1 - res_vec1.min(axis=0)) / (res_vec1.max(axis=0) - res_vec1.min(axis=0))) * (
				range[1] - range[0]) + range[0]
	res_vec2_scaled = ((res_vec2 - res_vec2.min(axis=0)) / (res_vec2.max(axis=0) - res_vec2.min(axis=0))) * (
				range[1] - range[0]) + range[0]

	# print performance evaluation
	print("\nbest predicted rank in true:\n" + "mean:", np.mean(res_vec1_scaled), ", median:",
		  np.median(res_vec1_scaled))
	print("\nbest true rank in pred :\n" + "mean:", np.mean(res_vec2_scaled), ", median:",
		  np.median(res_vec2_scaled))

	# how many trees we need to check on average?
	res_vec2_scaled.sort()
	sorted_true_scaled = res_vec2_scaled[::-1]
	index99 = int(0.01 * len(sorted_true_scaled) - 1)  # index for loc 0.01 = 99%
	index95 = int(0.05 * len(sorted_true_scaled) - 1)
	#print("\nIn 0.99: {}%".format(sorted_true_scaled[index99]))
	print("\nIn 0.95: {}%".format(sorted_true_scaled[index95]))


	return





if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='arrange data for learning and implement learning algo')
	parser.add_argument('--step_number', '-st', required=True) 	 # counting from 1
	parser.add_argument('--operator_compose', '-op', default="add")
	args = parser.parse_args()

	dirpath = SUMMARY_FILES_DIR if platform.system() == 'Linux' else DATA_PATH
	df_orig = pd.read_csv(dirpath + CHOSEN_DATASETS_FILENAME)
	#df_prune = pd.read_csv(dirpath + "with_preds_prune.csv").reset_index(drop=True)   # _rank
	#df_rgft = pd.read_csv(dirpath + "with_preds_rgft.csv")  .reset_index(drop=True)   # _rank
	df_prune = pd.read_csv(dirpath + "with_preds_prune.csv").reset_index(drop=True)  # _rank
	df_rgft = pd.read_csv(dirpath + "with_preds_all_moves_rgft.csv").reset_index(drop=True)

	df_prune, df_rgft = preds_to_probs(df_prune, df_rgft)
	df_comp = compose_learners(df_prune, df_rgft, operator=args.operator_compose)
	rank_pred_by_ds, rank_test_by_ds = performance_evaluation(df_comp)
	#df_comp.to_csv(dirpath + "preds_ranked.csv")
	print_results(rank_pred_by_ds, rank_test_by_ds)


	#df_orig = df_orig[df_orig["path"].isin(paths_lst)]
	#df_orig["rank_pred"] = rank_pred_by_ds
	#df_orig[""]