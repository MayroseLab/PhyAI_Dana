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
		stats_filpath_for_params =SEP.join([SEP.join(tree_file.split(SEP)[:-4]), PHYML_STATS_FILENAME.format('bionj')])
		params_dict = parse_phyml_stats_output(msa_file_full_path, stats_filpath_for_params)
		f = [params_dict["fA"],params_dict["fC"],params_dict["fG"],params_dict["fT"]]
		rates = [params_dict["subAC"], params_dict["subAG"], params_dict["subAT"], params_dict["subCG"], params_dict["subCT"], params_dict["subGT"]]

		exec_line = run_phyml_TRUEmodelparams(msa_file_full_path, tree_file, base_model, params_dict["gamma"], params_dict["pInv"], f, rates)
		return exec_line
	else:
		execution_tags = " ".join(PHYML_MODEL_TAGS[base_model] + [PHYML_PINV_TAGS[pinv], PHYML_GAMMA_TAGS[gamma],
							   PHYML_TOPOLOGY_TAGS[topology_tag], PHYML_GENERAL_TAGS, "--run_id " + run_id])
		if tree_file:
			execution_tags += " -u " + tree_file
		return " ".join([PHYML_SCRIPT, "-i", msa_file_full_path, execution_tags])


def create_phyml_exec_line_full_model(msa_file_full_path, full_model, topology_tag, tree_file=False):
	pinv = "+I" in full_model
	gamma = "+G" in full_model
	base_model = re.sub(r'\+.*', '', full_model)

	return create_phyml_exec_line(msa_file_full_path, base_model, pinv, gamma, topology_tag, tree_file)



def run_phyml_TRUEmodelparams(msa_filepath, tree_filepath, model_name, alpha, pinv, f, rates):
	'''
	running phyml in interactive mode to be able to fix substitution model parameters
	'''
	run_id = '\n'.join(['R','br'])    #
	subs_model = re.match("([^+]+).*", model_name).group(1)
	models_dict_fixed = {'JC': '000000', 'F81': '000000', 'K80': '010010', 'HKY': '010010', 'SYM': '012345','GTR': '012345'}
	model_value = '\n'.join(['M\nM\nM\nM\nK', str(models_dict_fixed.get(subs_model))])

	f_interactive = 'E\n' + "\n".join(f)
	rates_interactive = "\n".join(rates)
	alpha_interactive = 'R\n' if not alpha else 'C\n4\nA\nn\n' + alpha
	pinv_interactive = '' if not pinv else 'V\nn\n' + pinv
	if_replace_output = "\n".join(['R','R'])   # R for replace, A for append
	if_subs_opt = "O"        # O for turning off, None for leaving default ml optimisation
	if_topology_opt = "O"    # O for turning off, None for leaving default (NNI) ml optimisation
	input_tree_mode = "U"    # U for changing to user tree input, None for leaving default bionj starting tree
	if_LRT = '\n'.join(['A', 'A'])

	params_interactive = '\n'.join([str(msa_filepath), if_replace_output, run_id, "+", model_value, rates_interactive, f_interactive, alpha_interactive, pinv_interactive, if_subs_opt, "+", if_topology_opt, input_tree_mode, "+", if_LRT, "Y", tree_filepath, ""])
	#params_interactive = str(msa_filepath) + '\nA\nA\nR\n' + run_id + "\n+\nM\nM\nM\nM\nK\n" + str(model_value) + "\n" + rates_interactive + f_interactive + alpha_interactive + pinv_interactive + "O\n+\n+\nA\nA\nY\n"
	interactiveM_command = "echo '" + params_interactive + "' | " + PHYML_SCRIPT + "\n"

	os.system("#!/bin/tcsh\n\n")
	#os.system(interactiveM_command)

	return interactiveM_command



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

