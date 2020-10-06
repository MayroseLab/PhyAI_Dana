import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")
import warnings

from defs import *


df1 = pd.read_csv("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/summary_files/with_preds_merged_20_1_example_5000_ytransformed_exp.csv")
#df2 = pd.read_csv("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/summary_files/with_preds_merged_20_1_exampleV2_V2_ytransformed_exp.csv")
df2 = pd.read_csv("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/summary_files/with_preds_merged_20_1_exampleV2_V4_ytransformed_exp.csv")

df = df1.merge(df2,on="", left_index=True, right_index=True, suffixes=('_prune', '_rgft'))