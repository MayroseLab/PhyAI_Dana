import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")

from defs import *
from data_processing.traverse_data_dirs import traverse_data_dirs
from my_utils.msa_functions import *
from execute_programs.RAxML_NG import run_raxml



	

def generate_random_trees(datapath):
	msa_path = datapath + MSA_PHYLIP_FILENAME
	tips_names = get_seqs_dict(msa_path).keys()
	max_edge_length = random.random()  # a float between 0.1 and 1
	
	t = Tree()
	t.populate(len(tips_names), names_library=tips_names, random_branches=True, branch_range=((0, max_edge_length)))
	t.unroot()
	
	dest_path = datapath + RANDOM_TREE_DIRNAME
	dest_tree = dest_path + RANDOM_TREE_FILENAME
	if not os.path.exists(dest_path):
		os.mkdir(dest_path)
	t.write(format=1, outfile=dest_tree)
	shutil.copy(msa_path, dest_path + MSA_PHYLIP_FILENAME)
	
	phyml_for_ll(dest_path, opt='br', treepath=dest_tree)
	


def phyml_for_ll(dataset_path, opt="bionj", treepath=None):
	msa_filepath = dataset_path + MSA_PHYLIP_FILENAME
	job_name = "phyml_" + "_".join(dataset_path.split("/")[-2:-1])
	cmd = "python " + CODE_PATH + "execute_programs/Phyml.py " + "-f " + msa_filepath + " -br {}".format(opt)
	if treepath:
		cmd += " -t " + treepath
	#create_job_file.main(command=cmd, dirpath=dataset_path, sh_file=job_name + ".sh", multiply_jobs=False, priority=-1,
	#					 job_name=job_name)
	os.system(cmd)



def create_job_starting_tree(dataset_path, tree_type):
	print("**************************************\n", dataset_path)
	job_name = "{}_starting_tree_".format(tree_type) + "_".join(dataset_path.split("/")[-2:-1])
	cmd = "python " + CODE_PATH + "data_processing/create_starting_trees.py -ds " + dataset_path + " -ttype " + tree_type
	create_job_file.main(command=cmd, dirpath=dataset_path, sh_file=job_name + ".sh", multiply_jobs=False,
						 priority=-1, job_name=job_name)
	

def generate_starting_trees(datapath, tree_type):
	if not os.path.exists(datapath + MSA_PHYLIP_FILENAME):
		mask_species_names_in_msa(datapath + MSA_PHYLIP_FILENAME_NOT_MASKED, datapath + MSA_PHYLIP_FILENAME, datapath)

	if tree_type == 'random':
		generate_random_trees(datapath)
	if tree_type == 'bionj':
		phyml_for_ll(datapath, tree_type)
	if tree_type == 'parsimony':
		run_raxml(msa_path=datapath+MSA_PHYLIP_FILENAME, tree_path=None, mode='starting_optimized', runover=True)




if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='perform all SPR moves')
	parser.add_argument('--dataset_path', '-ds', default=None)
	parser.add_argument('--tree_type', '-ttype', default=None)  # could be bionj | random | parsimony
	parser.add_argument('--index_to_start_run', '-istart', default=False)
	parser.add_argument('--nline_to_run', '-nlines', default=False)
	args = parser.parse_args()

	datapath = args.dataset_path
	if datapath:
		generate_starting_trees(datapath, args.tree_type)
	else:
		csv_path = SUMMARY_FILES_DIR + CHOSEN_DATASETS_FILENAME
		traverse_data_dirs(create_job_starting_tree, csv_path, (args.index_to_start_run, args.nline_to_run), args.tree_type)