import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")


from defs import *

from subprocess import Popen, PIPE, STDOUT
RAXML_NG_SCRIPT = "raxml-ng"

GROUP_ID = 'group_id'



def rf_rank(df_by_ds, sortby, location):
	best_pred_ix = df_by_ds[sortby].idxmax()



	temp_df = df_by_ds.sort_values(by=location, ascending=False).reset_index()  # changed ascending to False
	best_pred_rank = min(temp_df.index[temp_df["index"] == best_pred_ix].tolist())
	best_pred_rank += 1  # convert from pythonic index to position


	return best_pred_rank


def index_rf_scores(df_by_ds, t1):
	trees_df_ds = pd.read_csv(TREES_PER_DS.format(df_by_ds['path'].values[0], '1'))

	for i, row in trees_df_ds.iterrows():
		tree_str = row['newick']
		t2 = Tree(tree_str, format=1)
		if ROOTLIKE_NAME + "_2" in t2.get_leaf_names():
			t2.delete(t2 & ROOTLIKE_NAME + "_2")
		else:
			t2.set_outgroup(t2 & 'Sp0000')

		rf_score = t1.robinson_foulds(t2)[0]
		trees_df_ds.loc[i, "rf"] = trees_df_ds

		print(trees_df_ds['rf'])
		exit()




if __name__ == '__main__':
	df = pd.read_csv("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/data_test/with_preds_test_subs.csv")
	ml_tree_str = "(Sp0019:0.01073154,Sp0023:0.00744136,(((Sp0025:0.16862802,(Sp0027:0.13785334,Sp0011:0.17851188):0.01217161):0.01254977,((Sp0010:0.12270874,(Sp0012:0.12414156,(Sp0004:0.11076039,Sp0005:0.09839716):0.08062449):0.05245152):0.02885575,((Sp0008:0.07917821,Sp0020:0.08648837):0.11337615,(Sp0002:0.18288849,Sp0001:0.12425424):0.02513772):0.02860382):0.01285754):0.01753573,(((Sp0007:0.09300627,Sp0016:0.04971814):0.06045486,(Sp0022:0.18892855,(Sp0003:0.15118672,(Sp0000:0.09667597,Sp0015:0.08336488):0.01684364):0.02471939):0.02186005):0.02225697,((Sp0013:0.29550040,(Sp0026:0.28248980,Sp0021:0.27621247):0.04592804):0.04467169,(Sp0014:0.14784391,(Sp0006:0.31906820,((Sp0017:0.21935344,Sp0018:0.18520245):0.10051093,(Sp0024:0.24090004,Sp0009:0.26059768):0.05867151):0.03309910):0.05972386):0.02810186):0.02420151):0.01896213):0.23098854);"
	ml_tree = Tree(ml_tree_str, format=1)
	ml_tree.set_outgroup(ml_tree&'Sp0000')

	rank_pred_by_ds, rank_test_by_ds = {}, {}
	label = LABEL.format('merged')
	grouped_df_by_ds = df.groupby(FEATURES[GROUP_ID], sort=False)
	for group_id, df_by_ds in grouped_df_by_ds:
		index_rf_scores(df_by_ds, ml_tree)

		rank_pred_by_ds[group_id] = rf_rank(df_by_ds, "pred", label)
		rank_test_by_ds[group_id] = rf_rank(df_by_ds, label, "pred")




	'''
	ml_tree = "(Sp0019:0.01073154,Sp0023:0.00744136,(((Sp0025:0.16862802,(Sp0027:0.13785334,Sp0011:0.17851188):0.01217161):0.01254977,((Sp0010:0.12270874,(Sp0012:0.12414156,(Sp0004:0.11076039,Sp0005:0.09839716):0.08062449):0.05245152):0.02885575,((Sp0008:0.07917821,Sp0020:0.08648837):0.11337615,(Sp0002:0.18288849,Sp0001:0.12425424):0.02513772):0.02860382):0.01285754):0.01753573,(((Sp0007:0.09300627,Sp0016:0.04971814):0.06045486,(Sp0022:0.18892855,(Sp0003:0.15118672,(Sp0000:0.09667597,Sp0015:0.08336488):0.01684364):0.02471939):0.02186005):0.02225697,((Sp0013:0.29550040,(Sp0026:0.28248980,Sp0021:0.27621247):0.04592804):0.04467169,(Sp0014:0.14784391,(Sp0006:0.31906820,((Sp0017:0.21935344,Sp0018:0.18520245):0.10051093,(Sp0024:0.24090004,Sp0009:0.26059768):0.05867151):0.03309910):0.05972386):0.02810186):0.02420151):0.01896213):0.23098854);"
	best_pred_tree = "(Sp0023:0.00744136,((((Sp0010:0.122709,(Sp0012:0.124142,Sp0005:0.179022)N21:0.0524515)N12:0.0288557,((Sp0027:0.137853,Sp0011:0.178512)N11:0.0121716,((Sp0008:0.0791782,Sp0020:0.0864884)N22:0.113376,(Sp0002:0.182888,Sp0001:0.124254)N23:0.0251377)N13:0.0143019)N6:0.0143019)N7:0.0128575,Sp0025:0.181178)N4:0.0175357,(((Sp0007:0.0930063,Sp0016:0.0497181)N14:0.0604549,(Sp0022:0.188929,(Sp0003:0.151187,(Sp0000:0.096676,Sp0015:0.0833649)N39:0.0168436)N27:0.0247194)N15:0.02186)N8:0.022257,((Sp0013:0.2955,(Sp0026:0.28249,Sp0021:0.276212)N29:0.045928)N16:0.0446717,(Sp0014:0.147844,(Sp0006:0.319068,((Sp0017:0.219353,Sp0018:0.185202)N48:0.100511,(Sp0024:0.2409,Sp0009:0.260598)N49:0.0586715)N43:0.0330991)N31:0.0597239)N17:0.0281019)N9:0.0242015)N5:0.0189621)N3:0.230989,(Sp0004:0.11076,Sp0019:0.00536575)N33:0.00536575);"

	t1 = Tree(ml_tree, format=1)
	t2 = Tree(best_pred_tree, format=1)
	t1.set_outgroup(t1&'Sp0000')
	if ROOTLIKE_NAME + "_2" in t2.get_leaf_names():
		t2.delete(t2&ROOTLIKE_NAME + "_2")
	else:
		t2.set_outgroup(t2&'Sp0000')

	res = t1.robinson_foulds(t2)
	rf_score = res[0]
	max_rf = res[1]
	print(rf_score, max_rf)
	'''
