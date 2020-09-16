import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")

from defs import *
from parsing.parse_phyml import parse_phyml_stats_output


###############################
RAXML_NG_SCRIPT = "raxml-ng"
### ## this line is not general - fixed for: bl otimization (only) | GTR+I+G model | fixed-extracted-model-params |  #####
RAXML_CMD = RAXML_NG_SCRIPT + " --evaluate --msa {msa_file} --opt-branches on --opt-model off --model GTR{rates}+I{pinv}+G{alpha}+F{freq} --tree {tree}"# --nofiles interim --log RESULT" # --nofiles"
###############################


def run_raxml(msa_path, tree_path, runover=False):
	rates, pinv, alpha, freq = extract_model_params(msa_path, tree_path)

	RAxML_cmd = RAXML_CMD.format(msa_file=msa_path, tree=tree_path,
								 rates="{{{0}}}".format("/".join(rates)),
								 pinv="{{{0}}}".format(pinv), alpha="{{{0}}}".format(alpha),
								 freq="{{{0}}}".format("/".join(freq)))
	if runover:
		RAxML_cmd += " --redo"

	cmd = RAxML_cmd + "\n"
	os.system(cmd)


def extract_model_params(msa_file_full_path, tree_file):
	#stats_filpath_for_params = SEP.join([SEP.join(tree_file.split(SEP)[:-4]), PHYML_STATS_FILENAME.format('bionj')])
	# todo : make more generic before running on all!
	stats_filpath_for_params = SEP.join([SEP.join(tree_file.split(SEP)[:-1]), PHYML_STATS_FILENAME.format('bionj')])
	params_dict = parse_phyml_stats_output(msa_file_full_path, stats_filpath_for_params)
	freq = [params_dict["fA"], params_dict["fC"], params_dict["fG"], params_dict["fT"]]
	rates = [params_dict["subAC"], params_dict["subAG"], params_dict["subAT"], params_dict["subCG"],
			 params_dict["subCT"], params_dict["subGT"]]

	return rates, params_dict["pInv"], params_dict["gamma"], freq




if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='run raxml-ng')
	parser.add_argument('--msa_filepath', '-f', default=None)
	parser.add_argument('--tree_filepath', '-t', default=None)
	parser.add_argument('--runover', '-r', default=False, action='store_true')
	args = parser.parse_args()

	import time
	start_time = time.time()
	run_raxml(args.msa_filepath, args.tree_filepath, args.runover)
	print(time.time() - start_time)


