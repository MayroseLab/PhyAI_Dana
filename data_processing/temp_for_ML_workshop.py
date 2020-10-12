import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")
import warnings
warnings.filterwarnings("ignore")			# TEMP

from defs import *
from utils.msa_functions import *
from data_processing.create_starting_trees import *
SUBMISSION_DIR = DIRPATH + "submission_data/"




def copy_datasets(df, set):
	for index, row in df.iterrows():
		dataset_path = row["path"]
		dataset_name = re.search("(\/.*data2?\/.*\/)(.*\/.*)", row["path"]).group(2)
		dest_dirname = dataset_path.split("/")[-2]
		
		orig_msa = dataset_path + SEP + MSA_PHYLIP_FILENAME
		#dest_path_dir = SEP.join([SUBMISSION_DIR + set, dest_dirname, ""])
		dest_path_dir = SEP.join(["/groups/itay_mayrose/danaazouri/PhyAI/ML_workshop/data/" + set, dest_dirname, ""])
		if not os.path.exists(dest_path_dir):
			os.makedirs(dest_path_dir)

		shutil.copy(orig_msa, dest_path_dir + MSA_PHYLIP_FILENAME)
		#if

		exit()




if __name__ == '__main__':
	#df_validation = pd.read_csv("/groups/itay_mayrose/danaazouri/PhyAI/validation_set2/summary_files/" + CHOSEN_DATASETS_FILENAME)
	#df_6000 = pd.read_csv(SUMMARY_FILES_DIR + CHOSEN_DATASETS_FILENAME)
	#df_3700 = pd.read_csv(SUMMARY_FILES_DIR + "scores_per_ds_28_1_ytransformed_exp_minus.csv")
	#df_additional = df_6000[~df_6000["path"].isin(df_3700["path"].unique())]

	df = pd.read_csv(SUMMARY_FILES_DIR + CHOSEN_DATASETS_FILENAME)
	copy_datasets(df, "training_datasets")

	#copy_datasets(df_additional, "additional_training_data")
	#copy_datasets(df_validation, "validation_data")
	
	






