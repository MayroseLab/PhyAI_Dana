import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")

from defs import *

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

	return rf_of_best_pred-df_by_ds_with_rf['rf'].min(), ml_rank_in_pred


def index_rf_scores(df_by_ds, t1):
	trees_df_ds = pd.read_csv(TREES_PER_DS.format(df_by_ds['path'].values[0], '1'))
	trees_df_ds = trees_df_ds[trees_df_ds['rgft_name'] != SUBTREE1]
	trees_df_ds = trees_df_ds[trees_df_ds['rgft_name'] != SUBTREE2]

	for i, row in trees_df_ds.iterrows():
		tree_str = row['newick']
		t2 = Tree(tree_str, format=1)
		if ROOTLIKE_NAME + "_2" in t2.get_leaf_names():
			t2.delete(t2 & ROOTLIKE_NAME + "_2")
		else:
			t2.set_outgroup(t2 & 'Sp0000')  # todo: change back to Sp0000 (4 zeros)

		rf_score = t2.robinson_foulds(t1)[0]
		df_by_ds.loc[(df_by_ds['prune_name'] == row['prune_name']) & (df_by_ds['rgft_name'] == row['rgft_name']), "rf"] = rf_score

	return df_by_ds





if __name__ == '__main__':
	df = pd.read_csv("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/summary_files/with_preds_merged_20_1_ml_minus1_set_4200_ytransformed_exp.csv")
	ml_tree_str = "(Sp0019:0.01073154,Sp0023:0.00744136,(((Sp0025:0.16862802,(Sp0027:0.13785334,Sp0011:0.17851188):0.01217161):0.01254977,((Sp0010:0.12270874,(Sp0012:0.12414156,(Sp0004:0.11076039,Sp0005:0.09839716):0.08062449):0.05245152):0.02885575,((Sp0008:0.07917821,Sp0020:0.08648837):0.11337615,(Sp0002:0.18288849,Sp0001:0.12425424):0.02513772):0.02860382):0.01285754):0.01753573,(((Sp0007:0.09300627,Sp0016:0.04971814):0.06045486,(Sp0022:0.18892855,(Sp0003:0.15118672,(Sp0000:0.09667597,Sp0015:0.08336488):0.01684364):0.02471939):0.02186005):0.02225697,((Sp0013:0.29550040,(Sp0026:0.28248980,Sp0021:0.27621247):0.04592804):0.04467169,(Sp0014:0.14784391,(Sp0006:0.31906820,((Sp0017:0.21935344,Sp0018:0.18520245):0.10051093,(Sp0024:0.24090004,Sp0009:0.26059768):0.05867151):0.03309910):0.05972386):0.02810186):0.02420151):0.01896213):0.23098854);"
	ml_tree = Tree(ml_tree_str, format=1)
	ml_tree.set_outgroup(ml_tree&'Sp0000')
	#df = pd.read_csv("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/data_test/with_preds_test_subs.csv")
	#ml_tree_str = "(((((((Sp012:0.0468574,((Sp010:0.0185394,Sp003:0.0383335)N22:0.00450131,Sp002:0.0407356)N19:0.0119676)N14:1.3e-07,((((Sp026:0.084224,Sp005:0.0905619)N30:0.0208792,Sp017:0.0994158)N24:0.00190321,Sp015:0.0833994)N20:0.0203127,((((((Sp025:0.00209632,(Sp024:0.00209495,Sp022:0.00851277)N49:8e-08)N44:0.00202704,Sp008:0.00650192)N40:0.00436954,Sp019:0.00434327)N38:0.00203008,(Sp014:0.012611,(Sp001:0.00623502,Sp000:0.00206603)N43:9e-08)N39:0.0244698)N32:0.012999,Sp016:0.0455011)N26:0.0421575,(Sp009:0.0176727,Sp007:0.0549442)N27:0.00833237)N21:0.00275245)N15:0.0303858)N10:1.9e-07,Sp023:0.0926691)N8:0.0147304,((Sp020:0.108462,Sp018:0.0713038)N12:0.0243208,Sp011:0.0510737)N9:1.11e-06)N6:0.0178948,Sp004:0.0699032)N4:0.00756702,Sp013:0.055996)N1:0.0364941,Sp021:0.0101297,Sp006:0.0228496);"
	#ml_tree = Tree(ml_tree_str, format=1)
	#ml_tree.set_outgroup(ml_tree & 'Sp000')

	rf_best_pred_by_ds, rank_rf0_by_ds = {}, {}
	label = LABEL.format('merged')
	grouped_df_by_ds = df.groupby(FEATURES[GROUP_ID], sort=False)
	for group_id, df_by_ds in grouped_df_by_ds:
		df_by_ds_with_rf = index_rf_scores(df_by_ds, ml_tree)

		rf_best_pred_by_ds[group_id], rank_rf0_by_ds[group_id] = rf_rank(df_by_ds_with_rf)
		print(rf_best_pred_by_ds[group_id], rank_rf0_by_ds[group_id])
		exit()
