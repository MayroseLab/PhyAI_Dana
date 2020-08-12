import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")

from utils.general_utils import change_path_permissions_to_777
from defs import *
from parsing.parse_phyml import parse_phyml_stats_output


############################### additional parameters ###############################
PHYML_PINV_TAGS = {True: "-v e",
				   False: ""}
PHYML_GAMMA_TAGS = {True: "-a e -c 4",
					False: "-c 1"}
PHYML_TOPOLOGY_TAGS = {"ml": "-o tlr -s NNI",
					   "bionj": "-o lr",
                       "br": "-o l",
                       "no_opt": "-o n"}

################################# model parameters #################################
PHYML_BASE_FREQS_TAGS = {True: "-f m",
						 False: "-f 0.25,0.25,0.25,0.25"}
PHYML_SUBS_RATES_TAGS = {1: "-m 000000",
						 2: "-m 010010",
						 3: "-m 012345"}
PHYML_MODEL_TAGS = {"JC": [PHYML_SUBS_RATES_TAGS[1], PHYML_BASE_FREQS_TAGS[False]],
					"F81": [PHYML_SUBS_RATES_TAGS[1], PHYML_BASE_FREQS_TAGS[True]],
					"K80": [PHYML_SUBS_RATES_TAGS[2], PHYML_BASE_FREQS_TAGS[False]],
					"HKY": [PHYML_SUBS_RATES_TAGS[2], PHYML_BASE_FREQS_TAGS[True]],
					"SYM": [PHYML_SUBS_RATES_TAGS[3], PHYML_BASE_FREQS_TAGS[False]],
					"GTR": [PHYML_SUBS_RATES_TAGS[3], PHYML_BASE_FREQS_TAGS[True]]}

###################################### general ######################################
PHYML_GENERAL_TAGS = "-d nt -n 1 -b 0 --no_memory_check"
PHYML_SCRIPT = "/groups/itay_mayrose/shiranabad/applications/jmodeltest-2.1.7/exe/phyml/PhyML_3.0_linux64"


def create_phyml_exec_line(msa_file_full_path, base_model, pinv, gamma, topology_tag, tree_file=False):
	run_id = topology_tag

	# enable fixed model parameters
	if topology_tag == "br":
		stats_filpath_for_params =SEP.join([SEP.join(tree_file.split(SEP)[:-1]), PHYML_STATS_FILENAME.format('bionj')])
		params_dict = parse_phyml_stats_output(msa_file_full_path, stats_filpath_for_params)
		print(params_dict["pInv"], params_dict["gamma"],params_dict["fA"],params_dict["fC"],params_dict["fG"],params_dict["fT"])
		pinv_tags = "-v {}".format(params_dict["pInv"])
		gamma_tags = "-a {} -c 4".format(params_dict["gamma"])
		freq_tags = "-f {},{},{},{}".format(params_dict["fA"],params_dict["fC"],params_dict["fG"],params_dict["fT"])   # "fA", "fC", "fG", "fT"
		model_tags = [PHYML_SUBS_RATES_TAGS[3], freq_tags]
	else:
		pinv_tags = PHYML_PINV_TAGS[pinv]
		gamma_tags = PHYML_GAMMA_TAGS[gamma]
		model_tags = PHYML_MODEL_TAGS[base_model]

	execution_tags = " ".join(model_tags + [pinv_tags, gamma_tags,
							   PHYML_TOPOLOGY_TAGS[topology_tag], PHYML_GENERAL_TAGS, "--run_id " + run_id])
	if tree_file:
		execution_tags += " -u " + tree_file

	return " ".join([PHYML_SCRIPT, "-i", msa_file_full_path, execution_tags])


def create_phyml_exec_line_full_model(msa_file_full_path, full_model, topology_tag, tree_file=False):
	pinv = "+I" in full_model
	gamma = "+G" in full_model
	base_model = re.sub(r'\+.*', '', full_model)

	return create_phyml_exec_line(msa_file_full_path, base_model, pinv, gamma, topology_tag, tree_file)


def run_phyml(msa_filepath, full_model, topology_tag, tree_file=None, force_run=False, copy_msa=False):
	if copy_msa:
		msa_copy_path = SEP.join([SEP.join(tree_file.split(SEP)[:-1]), MSA_PHYLIP_FILENAME])
		shutil.copy(msa_filepath, msa_copy_path)
		msa_filepath = msa_copy_path
	phyml_exec_line = create_phyml_exec_line_full_model(msa_filepath, full_model, topology_tag, tree_file=tree_file)
	output_filename = "{}_phyml_{}_" + topology_tag + ".txt"
	
	if not force_run and os.path.exists(output_filename.format(msa_filepath, "stats")):
			return output_filename.format(msa_filepath, "stats")
	os.system(phyml_exec_line)


	return output_filename.format(msa_filepath, "{}")






if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='run phyml')
	parser.add_argument('--msa_filepath', '-f', default=None)
	parser.add_argument('--tree_filepath', '-t', default=None)
	parser.add_argument('--topology_tag', '-br', default="bionj")
	parser.add_argument('--model', '-m', default=MODEL_DEFAULT)
	parser.add_argument('--runover', '-r', default=False, action='store_true')
	parser.add_argument('--cpmsa', '-cp', default=False, action='store_true')
	args = parser.parse_args()

	output_filename = run_phyml(msa_filepath=args.msa_filepath, full_model=args.model,topology_tag=args.topology_tag,
								tree_file=args.tree_filepath, force_run=args.runover, copy_msa=args.cpmsa)
	print("DONE phyml execution:", output_filename.format("tree"))

