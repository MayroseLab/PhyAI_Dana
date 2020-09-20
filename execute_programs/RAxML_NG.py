import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")

from defs import *
from parsing.parse_phyml import parse_phyml_stats_output
from utils.msa_functions import *


###############################
RAXML_NG_SCRIPT = "raxml-ng"
RAXML_FIXED_SUBS = RAXML_NG_SCRIPT + " --evaluate --msa {msa_file} --opt-branches on --opt-model off --model GTR{rates}+I{pinv}+G{alpha}+F{freq} --tree {tree}"# --nofiles interim --log RESULT" # --nofiles"
RAXML_STARTING = RAXML_NG_SCRIPT + " --start pars{npars} --msa {msa_file} --opt-branches on --opt-model on --model GTR+I+G"
RAXML_PARTITIONS = RAXML_NG_SCRIPT + ""
###############################

PAR_LENGTH = 50



def extract_model_params(msa_file_full_path, tree_file):
	stats_filpath_for_params = SEP.join([SEP.join(tree_file.split(SEP)[:-4]), PHYML_STATS_FILENAME.format('bionj')])  # [-4] cause each tree_file (when running all SPR) is in /rearrangements/prune_name/rgft_name/tree_filename
	params_dict = parse_phyml_stats_output(msa_file_full_path, stats_filpath_for_params)
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
		''' 
		GTR+I+G model (fixed-extracted-params) | bl otimization 
		'''
		rates, pinv, alpha, freq = extract_model_params(msa_path, tree_path)
		RAxML_cmd = RAXML_FIXED_SUBS.format(msa_file=msa_path, tree=tree_path,
									 rates="{{{0}}}".format("/".join(rates)),
									 pinv="{{{0}}}".format(pinv), alpha="{{{0}}}".format(alpha),
									 freq="{{{0}}}".format("/".join(freq)))
	elif mode == 'parsimony-optimized':
		from Bio.Phylo.Applications import RaxmlCommandline
		#RAxML_cmd = RAXML_STARTING.format(msa_file=msa_path, npars='{1}')#, nrand='{1}')
		try:
			raxml_cline = RaxmlCommandline(sequences=msa_path, model="GTR+I+G", parsimony=True)
			out, err = raxml_cline()
		except:
			print("xxx")

	elif mode == 'partitions':
		partitions_dirpath = SEP.join([SEP.join(tree_path.split(SEP)[:-1]), "partitions/"])
		n_partitions = partition_msa(msa_path, partitions_dirpath)
		RAxML_cmd = ''
		# todo: continue
	else:
		RAxML_cmd = 'complete_if_needed'

	#if runover:
	#	RAxML_cmd += " --redo"
	#os.system(RAxML_cmd + "\n")







if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='run raxml-ng')
	parser.add_argument('--msa_filepath', '-f', default=None)
	parser.add_argument('--tree_filepath', '-t', default=None)
	parser.add_argument('--opt_mode', '-mo', default='fixed_subs')  # could be: fixed_subs | parsimony-optimized | partitions
	parser.add_argument('--runover', '-r', default=False, action='store_true')
	args = parser.parse_args()

	#exit()
	import time
	start_time = time.time()
	run_raxml(args.msa_filepath, args.tree_filepath, args.opt_mode, args.runover)
	print(time.time() - start_time)


