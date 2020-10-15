import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")
import warnings
warnings.filterwarnings("ignore")			# TEMP

from defs import *
from utils.msa_functions import *
from data_processing.create_starting_trees import *
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

SUBMISSION_DIR = DIRPATH + "submission_data/"



def trucate_to_20(msa):
	alignment = AlignIO.read(msa, PHYLIP_FORMAT)
	new_align = alignment[0:20, :]
	trunc_msa = MultipleSeqAlignment(new_align)
	AlignIO.write(trunc_msa, msa, PHYLIP_FORMAT)
	rewrite_in_phylip(msa)


def copy_datasets(df, set):
	new_df = pd.DataFrame()
	for index, row in df.iterrows():
		dataset_path = row["path"]
		orig_msa = dataset_path + SEP + MSA_PHYLIP_FILENAME
		ntaxa, nchars = get_msa_properties(orig_msa)
		if ntaxa < 20:
			continue

		datapath = "/groups/itay_mayrose/danaazouri/PhyAI/ML_workshop/data/" + set
		dest_dirname = dataset_path.split("/")[-2]
		dest_path_dir = SEP.join([datapath, dest_dirname, ""])
		new_df.loc[index, "path"] = dest_path_dir
		if not os.path.exists(dest_path_dir):
			os.makedirs(dest_path_dir)
		msa = dest_path_dir + MSA_PHYLIP_FILENAME
		shutil.copy(orig_msa, msa)

		if ntaxa > 20:
			trucate_to_20(msa)

		generate_starting_trees(dest_path_dir, 'bionj')

	new_df.to_csv("/groups/itay_mayrose/danaazouri/PhyAI/ML_workshop/datasets_20taxa.csv")






if __name__ == '__main__':
	#df_validation = pd.read_csv("/groups/itay_mayrose/danaazouri/PhyAI/validation_set2/summary_files/" + CHOSEN_DATASETS_FILENAME)
	#df_6000 = pd.read_csv(SUMMARY_FILES_DIR + CHOSEN_DATASETS_FILENAME)
	#df_3700 = pd.read_csv(SUMMARY_FILES_DIR + "scores_per_ds_28_1_ytransformed_exp_minus.csv")
	#df_additional = df_6000[~df_6000["path"].isin(df_3700["path"].unique())]

	df = pd.read_csv(SUMMARY_FILES_DIR + CHOSEN_DATASETS_FILENAME)
	copy_datasets(df, "training_datasets")

	#copy_datasets(df_additional, "additional_training_data")
	#copy_datasets(df_validation, "validation_data")
	
	






