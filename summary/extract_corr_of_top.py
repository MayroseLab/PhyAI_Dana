import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")

from defs import *

from statistics import mean, median
import warnings
warnings.filterwarnings("ignore")  # TEMP
pd.set_option('display.max_columns', 40)

GROUP_ID = 'group_id'




def ds_scores(df, p, move_type='merged'):
    label = LABEL.format(move_type)
    sp_corrs = []
    grouped_df_by_ds = df.groupby(FEATURES[GROUP_ID], sort=False)
    for i, (group_id, df_by_ds) in enumerate(grouped_df_by_ds):
        temp_df = df_by_ds.sort_values(by=label, ascending=False).reset_index().head(int(len(df_by_ds)*p))
        temp_df = temp_df[[label, "pred"]]

        sp_corr = temp_df.corr(method='spearman').iloc[1, 0]
        if sp_corr:
            #if sp_corr < 0:
            #    print(group_id, sp_corr)
            sp_corrs.append(sp_corr)
        else:
            print("XXXXX", group_id)
            sp_corrs.append(None)
    print("*****************************" , i)
    return sp_corrs



def print_and_index_results(df_datasets, sp_corrs_lst):
    #### score 1 ####
    spearman_corrs = sp_corrs_lst
    df_datasets['corr_p'] = spearman_corrs
    #print("\nsapearman corr top xx%:\n" + "mean:", mean([e for e in spearman_corrs if not math.isnan(e)]), ", median:",median(spearman_corrs))
    print("##########################")

    return df_datasets




if __name__ == '__main__':
    dirpath = SUMMARY_FILES_DIR if platform.system() == 'Linux' else DATA_PATH
    #df_with_preds = pd.read_csv("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/summary_files/with_preds_merged_19_1_4200_ytransformed_exp_1cpu.csv")
    df_with_preds = pd.read_csv("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/summary_files/results_saturation/with_preds_merged_19_1_4500_ytransformed_exp.csv")
    #df_datasets = pd.read_csv("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/summary_files/scores_per_ds_19_1_4200_ytransformed_exp_1cpu.csv")
    df_datasets = pd.read_csv("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/summary_files/results_saturation/scores_per_ds_19_1_4500_ytransformed_exp.csv")

    for p in [0.1, 0.25, 0.5, 0.75, 1]:
        sp_corrs_lst = ds_scores(df_with_preds, p)
        df_datasets = print_and_index_results(df_datasets, sp_corrs_lst)
        df_datasets.to_csv("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/summary_files/scores_per_ds_19_1_4200_top_{}%.csv".format(p*100))
        #print_and_index_results(df_datasets, sp_corrs_lst)

