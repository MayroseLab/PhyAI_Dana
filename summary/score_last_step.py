import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")

from defs import *
from utils.create_job_file import get_job_qsub_command

GROUP_ID = 'group_id'




def rf_rank(df_by_ds_with_rf):
	assert df_by_ds_with_rf['rf'].min() == 0.0
	best_pred_ix = df_by_ds_with_rf["pred"].idxmax()
	rf_of_best_pred = df_by_ds_with_rf.loc[best_pred_ix, 'rf']


	best_true_id = df_by_ds_with_rf[df_by_ds_with_rf['rf'] == 0.0]["pred"].idxmax()
	temp_df = df_by_ds_with_rf.sort_values(by='d_ll_merged', ascending=False).reset_index()
	ml_rank_in_pred = min(temp_df.index[temp_df["index"] == best_true_id].tolist())

	ml_rank_in_pred += 1  # convert from pythonic index to position
	ml_rank_in_pred /= len(df_by_ds_with_rf['rf'].index)  # scale the rank according to rankmax
	ml_rank_in_pred *= 100

	return rf_of_best_pred-(df_by_ds_with_rf['rf'].min()), ml_rank_in_pred


def index_rf_scores(df, path, t1):
	trees_df_ds = pd.read_csv(TREES_PER_DS.format(path, '1'))
	trees_df_ds = trees_df_ds[trees_df_ds['rgft_name'] != SUBTREE1]
	trees_df_ds = trees_df_ds[trees_df_ds['rgft_name'] != SUBTREE2]

	for i, row in trees_df_ds.iterrows():
		tree_str = row['newick']
		t2 = Tree(tree_str, format=1)
		if ROOTLIKE_NAME + "_2" in t2.get_leaf_names():
			t2.delete(t2 & ROOTLIKE_NAME + "_2")
		else:
			t2.set_outgroup(t2 & 'Sp0000')

		rf_score = t2.robinson_foulds(t1)[0]
		#df_by_ds.loc[(df_by_ds['prune_name'] == row['prune_name']) & (df_by_ds['rgft_name'] == row['rgft_name']), "rf"] = rf_score
		df.loc[(df['path'] == path) & (df['prune_name'] == row['prune_name']) & (df['rgft_name'] == row['rgft_name']), "rf"] = rf_score

	return df




def submit_job(id, start):
	print("**************************************   ", id, start)
	job_name = "calc_RF.sh"
	cmd = "python " + CODE_PATH + "summary/score_last_step.py -i " + start + " -id " + id

	qsub_cmd = get_job_qsub_command(job_name=job_name,
									command=cmd,
									error_files_path=SUMMARY_FILES_DIR + "error_files/")
	os.system(qsub_cmd)



if __name__ == '__main__':
	#'''
	parser = argparse.ArgumentParser(description='arrange data for learning and implement learning algo')
	parser.add_argument('--istart', '-i', default=None)
	parser.add_argument('--subset_id', '-id', default=None)
	args = parser.parse_args()
	
	size = 61880
	if not args.istart:
		for id,start in enumerate(range(0, 6187910, size)):
			submit_job(str(id), str(start))

	else:
		df = pd.read_csv("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/summary_files/with_preds_merged_20_1_ml_minus1_set_4200_ytransformed_exp.csv")
		withRFscore = "/groups/itay_mayrose/danaazouri/PhyAI/DBset2/summary_files/with_preds_merged_20_1_ml_minus1_set_with_RFscore.csv"
		ml_tree_str = "(Sp0019:0.01073154,Sp0023:0.00744136,(((Sp0025:0.16862802,(Sp0027:0.13785334,Sp0011:0.17851188):0.01217161):0.01254977,((Sp0010:0.12270874,(Sp0012:0.12414156,(Sp0004:0.11076039,Sp0005:0.09839716):0.08062449):0.05245152):0.02885575,((Sp0008:0.07917821,Sp0020:0.08648837):0.11337615,(Sp0002:0.18288849,Sp0001:0.12425424):0.02513772):0.02860382):0.01285754):0.01753573,(((Sp0007:0.09300627,Sp0016:0.04971814):0.06045486,(Sp0022:0.18892855,(Sp0003:0.15118672,(Sp0000:0.09667597,Sp0015:0.08336488):0.01684364):0.02471939):0.02186005):0.02225697,((Sp0013:0.29550040,(Sp0026:0.28248980,Sp0021:0.27621247):0.04592804):0.04467169,(Sp0014:0.14784391,(Sp0006:0.31906820,((Sp0017:0.21935344,Sp0018:0.18520245):0.10051093,(Sp0024:0.24090004,Sp0009:0.26059768):0.05867151):0.03309910):0.05972386):0.02810186):0.02420151):0.01896213):0.23098854);"
		ml_tree = Tree(ml_tree_str, format=1)
		ml_tree.set_outgroup(ml_tree & 'Sp0000')

		id = args.subset_id
		start = int(args.istart)

		df_subs = df.iloc[start:start+int(size)]
		grouped_df_by_ds = df_subs.groupby(FEATURES[GROUP_ID], sort=False)
		for group_id, df_by_ds in grouped_df_by_ds:
			path = df_by_ds['path'].values[0]
			df_subs = index_rf_scores(df_subs, path, ml_tree)
			df_subs.to_csv("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/summary_files/with_preds_merged_20_1_ml_minus1_set_with_RFscore_subs{}.csv".format(id))
	'''

	df_path = "/groups/itay_mayrose/danaazouri/PhyAI/DBset2/summary_files/with_preds_merged_20_1_ml_minus1_set_with_RFscore.csv"
	if not os.path.exists(df_path):
		df = pd.DataFrame()
		for id in range(0, 100):
			subs_withRFscore = pd.read_csv("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/summary_files/with_preds_merged_20_1_ml_minus1_set_with_RFscore_subs{}.csv".format(id))
			df = pd.concat([df, subs_withRFscore], ignore_index=True)
		df.to_csv(df_path)
	else:
		df = pd.read_csv(df_path)

	#####################

	rf_best_pred_by_ds, rank_rf0_by_ds = {}, {}

	grouped_df_by_ds = df.groupby(FEATURES[GROUP_ID], sort=False)
	for group_id, df_by_ds in grouped_df_by_ds:
		path = df_by_ds['path'].values[0]
		rf_best_pred_by_ds[group_id], rank_rf0_by_ds[group_id] = rf_rank(df_by_ds)
		print(rf_best_pred_by_ds[group_id], rank_rf0_by_ds[group_id])

	from statistics import mean, median
	print(rf_best_pred_by_ds)
	print(rank_rf0_by_ds)
	print(mean(rf_best_pred_by_ds), median(rf_best_pred_by_ds))
	print(mean(rank_rf0_by_ds), median(rank_rf0_by_ds))
	'''
