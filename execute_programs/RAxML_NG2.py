import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")

from defs import *
from parsing.parse_phyml import parse_phyml_stats_output
from parsing.parse_raxml_NG import parse_raxmlNG_output
from utils.msa_functions import *


###############################
RAXML_NG_SCRIPT = "raxml-ng"
RAXML_EVALUATE = RAXML_NG_SCRIPT + " --evaluate --msa {msa_file} --threads 1 --opt-branches {opt_bl} --opt-model {opt_model} --model GTR{rates}+I{pinv}+G{alpha}+F{freq} --tree {tree}"# --nofiles interim --log RESULT" # --nofiles"
RAXML_STARTING = RAXML_NG_SCRIPT + " --tree {tree_type}{{1}} --msa {msa_file} --threads 1 --opt-branches off --opt-model off --model GTR+I+G --start"
RAXML_PARTITIONS = RAXML_NG_SCRIPT + ""
###############################

PAR_LENGTH = 50
SOFTWARE_STARTING_TREES = 'raxml'     # could be either 'raxml' or 'phyml
STARTING_TREE_TYPE='pars'             # could be either 'pars' or 'rand'


def extract_model_params(msa_file_full_path, tree_file, software):
	if software == 'phyml':
		stats_filpath_for_params = SEP.join([SEP.join(tree_file.split(SEP)[:-4]), PHYML_STATS_FILENAME.format('bionj')])  # [-4] cause each tree_file (when running all SPR) is in /rearrangements/prune_name/rgft_name/tree_filename
		params_dict = parse_phyml_stats_output(msa_file_full_path, stats_filpath_for_params)
	if software == 'raxml':
		log_filpath_for_params = SEP.join([SEP.join(tree_file.split(SEP)[:-4]), RAXML_STATS_FILENAME])
		#log_filpath_for_params = SEP.join([SEP.join(tree_file.split(SEP)[:-1]), 'real_msa4.phy.raxml.log'])
		params_dict = parse_raxmlNG_output(log_filpath_for_params)

	freq = [params_dict["fA"], params_dict["fC"], params_dict["fG"], params_dict["fT"]]
	rates = [params_dict["subAC"], params_dict["subAG"], params_dict["subAT"], params_dict["subCG"],
			 params_dict["subCT"], params_dict["subGT"]]

	return rates, params_dict["pInv"], params_dict["gamma"], freq


def partition_msa(msa_file, partitions_dirpath):
	"""
	splits the msa to chunks of 50's and saves each as a separate file
	"""
	# todo: just an outline. fit pipline<->script
	msa = get_msa_from_file(msa_file)
	ntaxa, nchars = get_msa_properties(msa)

	if not os.path.exists(partitions_dirpath + "part_%d.phy" % (nchars//PAR_LENGTH-1)):
		for i in range(0, nchars, PAR_LENGTH):
			cur_msa = msa[:, i:i+PAR_LENGTH]
			with open(partitions_dirpath + "part_%d.phy" % (i//PAR_LENGTH), "w") as fpw:
				AlignIO.write(cur_msa, fpw, PHYLIP_FORMAT)
	return (nchars//PAR_LENGTH)


def run_raxml(msa_path, tree_path, mode='fixed_subs', runover=False):

	if mode == 'fixed_subs':

		rates, pinv, alpha, freq = extract_model_params(msa_path, tree_path, software=SOFTWARE_STARTING_TREES)
		RAxML_cmd = RAXML_EVALUATE.format(msa_file=msa_path, tree=tree_path,
										  opt_bl='on', opt_model='off',
									 rates="{{{0}}}".format("/".join(rates)),
									 pinv="{{{0}}}".format(pinv), alpha="{{{0}}}".format(alpha),
									 freq="{{{0}}}".format("/".join(freq)))

	if mode == 'starting_optimized':
		## first generate a non-optimized parsimony tree
		RAxML_cmd_generate_start = RAXML_STARTING.format(msa_file=msa_path, tree_type=STARTING_TREE_TYPE)
		res = os.system(RAxML_cmd_generate_start + "\n")

		## then optimize all parameters
		starting_tree_path = msa_path + RAXML_STARTTREE_SUF
		if os.path.exists(starting_tree_path):
			RAxML_cmd = RAXML_EVALUATE.format(msa_file=msa_path, tree=starting_tree_path,
											  opt_bl='on', opt_model='on',
											  rates='', pinv='', alpha='', freq='')
		else:
			print(starting_tree_path, "does not exist\nthus starting tree optimization did NOT obtained")
			RAxML_cmd = ''

	if mode == 'partitions':
		partitions_dirpath = SEP.join([SEP.join(tree_path.split(SEP)[:-1]), "partitions/"])
		n_partitions = partition_msa(msa_path, partitions_dirpath)
		RAxML_cmd = ''
		# todo: continue


	if runover:
		RAxML_cmd += " --redo"
	res = os.system(RAxML_cmd + "\n")

	for suf in ['rba', 'bestModel', 'startTree', 'reduced.phy']:
		filepath = msa_path + '.raxml.' + suf
		if os.path.exists(filepath):
			os.remove(filepath)

	return






if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='run raxml-ng')
	parser.add_argument('--msa_filepath', '-f', default=None)
	parser.add_argument('--tree_filepath', '-t', default=None)
	parser.add_argument('--opt_mode', '-br', default='fixed_subs')  # could be: fixed_subs | starting_optimized | partitions
	parser.add_argument('--runover', '-r', default=False, action='store_true')
	parser.add_argument('--cpmsa', '-cp', default=False, action='store_true')
	args = parser.parse_args()

	if args.cpmsa:
		msa_copy_path = SEP.join([SEP.join(args.tree_filepath.split(SEP)[:-1]), MSA_PHYLIP_FILENAME])
		shutil.copy(args.msa_filepath, msa_copy_path)
		msa_filepath = msa_copy_path

	import time
	start_time = time.time()
	res = run_raxml(args.msa_filepath, args.tree_filepath, args.opt_mode, args.runover)
	print(time.time() - start_time)
