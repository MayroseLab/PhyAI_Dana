import sys, os
sys.path.append(os.path.dirname(sys.path[0]))

from defs import *
#from Bio import AlignIO
#from data_processing.create_starting_trees import phyml_for_ll

pd.set_option('display.max_columns', 10)




def copy_datasets(csv_file):
	df = pd.read_csv(csv_file)
	df_exist1 = pd.read_csv(dirpath + "sampled_datasets_Shiran_paths.csv")
	df_exist2 = pd.read_csv(dirpath + "sampled_datasets_Shiran_paths2.csv")
	
	df_exist1 = df_exist1[(df_exist1["Database"] == "orthoMam") | (df_exist1["Database"] == "PANDIT")]
	df = df[(df["Database"] == "orthoMam") | (df["Database"] == "PANDIT")]
	df = df[~df["path"].isin(df_exist1["path"].values)]
	df = df[~df["path"].isin(df_exist2["path"].values)]
	
	df = df.sample(n=155, random_state=1)
	df.to_csv(dirpath + "sampled_datasets_Shiran_paths3.csv")
	
	'''
	for index, row in df_exist.iterrows():
		dataset_name = re.search("(\/.*data2?\/.*\/)(.*\/.*)", row["path"]).group(2)
		name_split = dataset_name.split("/")
		dest_path_dir = DATA_PATH + name_split[0] + "_" + name_split[1] + SEP
		
		df_exist.loc[index, "path"] = dest_path_dir
		my_orig_path = SEP.join([DIRPATH + "validation_set", "data", "training_datasets", ""]) + name_split[0] + "_" + name_split[1] + SEP
		shutil.copytree(my_orig_path, dest_path_dir)
	'''
	for index, row in df.iterrows():
		#dataset_name = re.search("(\/.*)*((\/.*){2}\/)", row["path"]).group(2)
		dataset_name = re.search("(\/.*data2?\/.*\/)(.*\/.*)", row["path"]).group(2)
		name_split = dataset_name.split("/")
		dest_path_dir = DATA_PATH + name_split[0] + "_" + name_split[1] + SEP
		if not os.path.exists(dest_path_dir):
			os.mkdir(dest_path_dir)
		
		df.loc[index, "path"] = dest_path_dir
		shutil.copy(row["path"] + SEP + MSA_PHYLIP_FILENAME_NOT_MASKED, dest_path_dir + MSA_PHYLIP_FILENAME)
	
		phyml_for_ll(dest_path_dir)   # creates a job

	df = pd.concat([pd.read_csv(dirpath+CHOSEN_DATASETS_FILENAME), df])
	df.to_csv(dirpath + "sampled_datasets.csv")
	
	
def copy_gapped_datasets(dirpath):
	df_copy_from = pd.read_csv(dirpath + "gappad_datasets.csv")

	for index, row in df_copy_from.iterrows():
		# dataset_name = re.search("(\/.*)*((\/.*){2}\/)", row["path"]).group(2)
		dataset_name = re.search("(\/.*data2?\/.*\/)(.*\/.*)", row["path"]).group(2)
		name_split = dataset_name.split("/")
		dest_path_dir = DATA_PATH + name_split[0] + SEP
		if not os.path.exists(dest_path_dir):
			os.mkdir(dest_path_dir)

		df_copy_from.loc[index, "path"] = dest_path_dir
		shutil.copy(row["path"] + SEP + MSA_PHYLIP_FILENAME_NOT_MASKED, dest_path_dir + MSA_PHYLIP_FILENAME)

	df = pd.concat([pd.read_csv(dirpath + CHOSEN_DATASETS_FILENAME), df_copy_from], ignore_index=True)
	df.to_csv(dirpath + "sampled_datasets_updated.csv")




if __name__ == '__main__':
	dirpath = SUMMARY_FILES_DIR if platform.system() == 'Linux' else DATA_PATH
	#outpath = dirpath + CHOSEN_DATASETS_FILENAME
	#copy_datasets(dirpath + "validation_set_Shirnas_paths_7to70taxa.csv")

	copy_gapped_datasets(SUMMARY_FILES_DIR)
	
	'''
	df_exist = pd.read_csv(SEP.join([DIRPATH + "validation_set", "summary_files", ""]) + CHOSEN_DATASETS_FILENAME)
	df_exist = df_exist[(df_exist["Database"] != "orthoMam") & (df_exist["Database"] != "PANDIT") & (df_exist["Database"] != "Rfam")]
	for index, row in df_exist.iterrows():
		my_orig_path = row["path"]
		dest_path_dir = my_orig_path.replace("/validation_set/", "/DBset2/")
		df_exist.loc[index, "path"] = dest_path_dir
		shutil.copytree(my_orig_path, dest_path_dir)
	'''