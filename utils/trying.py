import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")

from defs import *








if __name__ == '__main__':
	df = pd.read_csv("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/summary_files/with_preds_merged_20_1_ml_minus1_set_4200_ytransformed_exp.csv")
	#best_pred_ix = df["pred"].idxmax()
	temp_df = df.sort_values(by='pred', ascending=False).reset_index()
	#print(temp_df['d_ll_merged'].head(20))
	print(temp_df['d_ll_prune'].head(30), "\n")

	temp_df = df.sort_values(by='d_ll_merged', ascending=False).reset_index()
	#print(temp_df['d_ll_merged'].head(20))
	print(temp_df['d_ll_prune'].head(30), "\n")
