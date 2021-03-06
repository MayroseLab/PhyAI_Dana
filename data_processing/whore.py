import sys, os
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")

from defs import *
from my_utils.tree_functions import *
from data_processing.traverse_data_dirs import traverse_data_dirs
from execute_programs.SPR_move import *

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

SUMMARIES_PER_DS_LST = ["ds_summary_prune_br_step1.csv", "ds_summary_prune_br_step2.csv", "ds_summary_rgft_br_step1.csv", "ds_summary_rgft_br_step2.csv", "newicks_step1.csv", "newicks_step2.csv"]
#SUMMARIES_PER_DS_LST_TEMP = ["newicks_step1.csv"] #, "ds_summary_rgft_br_step1.csv", "newicks_step1.csv"]
#SUMMARIES_PER_DS_LST_TEMP = ["newicks_step1.csv", "ds_summary_prune_br_step1.csv", "ds_summary_rgft_br_step1.csv", "masked_species_real_msa.phy_phyml_stats_bionj.txt", "masked_species_real_msa.phy_phyml_tree_bionj.txt", "masked_species_real_msa.phy_phyml_tree_bionj_no_internal.txt"]
SUMMARIES_PER_DS_LST_TEMP = ["ds_summary_rgft_br_step1.csv"]



def rearrange_dirs_for_rerun(datapath):
	new_dir = datapath + "v1/"
	if not os.path.exists(new_dir):
		os.mkdir(new_dir)
	os.system("mv " + datapath + "*.csv " + new_dir)
	#os.system("mv " + datapath + RAXML_TREE_FILENAME + "* " + new_dir)
	#os.system("mv " + datapath + RAXML_STATS_FILENAME + " " + new_dir)
	os.system("mv " + datapath + "*.txt " + new_dir)

	return


def copy_datasets(df_set, main_dir):
	for i, row in df_set.iterrows():
		path = row["path"]
		msa_orig = path + MSA_PHYLIP_FILENAME
		dest_dir = main_dir + path.split("/")[-2] + "/"

		if os.path.exists(msa_orig):
			if not os.path.exists(dest_dir):
				os.mkdir(dest_dir)
				shutil.copyfile(msa_orig, dest_dir + MSA_PHYLIP_FILENAME_NOT_MASKED)
		else:
			print(msa_orig, "does not exist !")
	return


def arrange_submission_data():
	df_all = pd.read_csv(SUMMARY_FILES_DIR + CHOSEN_DATASETS_FILENAME)
	df_4200 = pd.read_csv("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/summary_files/scores_per_ds_19_1_4200_ytransformed_exp_1cpu.csv")
	df_additional = df_all[~df_all["path"].isin(df_4200["path"].values)]
	dirpath = "/groups/itay_mayrose/danaazouri/PhyAI/submission_data/"

	copy_datasets(df_4200, dirpath + 'training_data/')
	copy_datasets(df_additional, dirpath + 'additional_training_data/')

	return


def delete_err_dirpath(datapath):
	err_dirpath = datapath + "error_files/"
	shutil.rmtree(err_dirpath, ignore_errors=True)
	#shutil.rmtree(datapath+REARRANGEMENTS_NAME+"s/", ignore_errors=True)
	#shutil.rmtree(datapath + RANDOM_TREE_DIRNAME, ignore_errors=True)


def make_sure_all_exist(datapath, missing_paths_per_ds):
	for filename in SUMMARIES_PER_DS_LST_TEMP:
		if not os.path.exists(datapath + filename):
			missing_paths_per_ds.append(datapath + filename)
	return missing_paths_per_ds


def missing_results():
	missing_paths_lst = []
	df = pd.read_csv(SUMMARY_FILES_DIR + CHOSEN_DATASETS_FILENAME)
	for index, row in df.iterrows():
		missing_paths_lst.extend(make_sure_all_exist(row["path"], []))

	with open(SUMMARY_FILES_DIR + "ds_to_rerun.txt", 'w') as fp:
		fp.write("\n".join(missing_paths_lst))
	return


def add_atts():
	df_scores = pd.read_csv(SUMMARY_FILES_DIR + SCORES_PER_DS.format("20_1_4300_ytransformed_exp"))
	tree_relpath = MSA_PHYLIP_FILENAME + "_phyml_tree_{0}.txt".format('bionj')
	for i,row in df_scores.iterrows():
		dirpath = row["path"]
		full_path = dirpath + tree_relpath
		tbl = get_total_branch_lengths(full_path)
		df_scores.ix[i, "tbl"] = tbl

	df_scores.to_csv(SUMMARY_FILES_DIR + SCORES_PER_DS.format("20_1_4300_ytransformed_exp_with_atts"))
	return

def create_data_dirs():
	df_paths = pd.DataFrame()
	df_trees = pd.read_csv("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/data/training_datasets/exampleSml/newicks_step1.csv")
	df_trees = df_trees[df_trees["rgft_name"] != SUBTREE1]
	df_trees = df_trees[df_trees["rgft_name"] != SUBTREE2]
	for i,row in df_trees.iterrows():
		tree_str = row['newick']
		dataset_dirpath = DATA_PATH + row['prune_name'] + '_' + row['rgft_name'] + SEP
		if not os.path.exists(dataset_dirpath):
			os.mkdir(dataset_dirpath)
		with open(dataset_dirpath + 'masked_species_real_msa.phy_phyml_tree_bionj.txt', 'w') as fp:
			fp.write(tree_str)

		df_paths.loc[i, "path"] = dataset_dirpath

	df_paths.to_csv(SUMMARY_FILES_DIR + CHOSEN_DATASETS_FILENAME)
	return


def do_something(datapath):
	#add_atts()
	#delete_err_dirpath(datapath)
	#rearrange_dirs_for_rerun(datapath)
	#missing_results()
	#create_data_dirs()
	arrange_submission_data()

	'''
	df = pd.read_csv(SUMMARY_FILES_DIR + LEARNING_DATA.format("all_moves", "1"))
	df = df.rename(columns={"res_tree_tbl": "res_tree_tbl_rgft", "res_tree_edge_length": "res_tree_edge_length_rgft"})
	print(df.columns)
	df.to_csv(SUMMARY_FILES_DIR + LEARNING_DATA.format("all_moves_up", "1"))
	'''




def create_job(dataset_path):
	job_name = "whore_job"
	cmd = "python " + CODE_PATH + "data_processing/whore.py -ds " + dataset_path
	res = create_job_file.main(command=cmd, dirpath=dataset_path, sh_file=job_name + ".sh", multiply_jobs=False,
						 priority=-1, job_name=job_name)


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--dataset_path', '-ds', default=None)
	parser.add_argument('--index_to_start_run', '-istart', default=False)
	parser.add_argument('--nline_to_run', '-nlines', default=False)
	args = parser.parse_args()

	datapath = args.dataset_path
	if datapath:
		do_something(datapath)
	else:
		csv_path = SUMMARY_FILES_DIR + CHOSEN_DATASETS_FILENAME
		traverse_data_dirs(create_job, csv_path, (args.index_to_start_run, args.nline_to_run))

