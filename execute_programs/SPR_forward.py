import sys

sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")
from execute_programs.SPR_move import *






def get_raxml_optimized_bl(tree_str, msa_file, rates, pinv, alpha, freq):
	model_line_params = 'GTR{rates}+I{pinv}+G{alpha}+F{freq}'.format(rates="{{{0}}}".format("/".join(rates)),
									 pinv="{{{0}}}".format(pinv), alpha="{{{0}}}".format(alpha),
									 freq="{{{0}}}".format("/".join(freq)))

	tree_rampath = "/dev/shm/" + str(random.random())  + str(random.random()) + "tree"  # the var is the str: tmp{dir_suffix}
	try:
		with open(tree_rampath, "w") as fpw:
			fpw.write(tree_str)

		p = Popen([RAXML_NG_SCRIPT, '--evaluate', '--msa', msa_file,'--threads', '2', '--opt-branches', 'on', '--opt-model', 'off', '--model', model_line_params, '--redo', '--tree', tree_rampath], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
		raxml_stdout = p.communicate()[0]
		raxml_output = raxml_stdout.decode()

		# parse to return also tree str
		ll = re.search("Final LogLikelihood:\s+(.*)", raxml_output).group(1).strip()
		with open(msa_file + ".raxml.bestTree", "r") as fpw:
			tree_str = fpw.read()


	except Exception as e:
		print(msa_file.split(SEP)[-1][3:])
		print(e)
		exit()
	finally:
		os.remove(tree_rampath)

	return ll, tree_str






if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='perform all SPR moves')
	parser.add_argument('--dataset_path', '-ds', default=None)
	parser.add_argument('--tree_relpath', '-trp', required=True)	# e.g., best_pred_st13
	args = parser.parse_args()

	dataset_path = args.dataset_path
	relpath = args.tree_relpath
	outpath = SUMMARY_PER_DS.format(dataset_path, "{}", "br", relpath)
	'''
	with open(dataset_path + relpath + '.txt', 'r') as fp:
		tree_str_not_opt = fp.read().strip()

	# update to string after bl optimization
	stats_filepath = dataset_path + PHYML_STATS_FILENAME.format("bionj")
	params_dict = (parse_phyml_stats_output(None,stats_filepath))
	freq, rates, pinv, alpha = [params_dict["fA"], params_dict["fC"], params_dict["fG"], params_dict["fT"]], [params_dict["subAC"], params_dict["subAG"], params_dict["subAT"], params_dict["subCG"], params_dict["subCT"],params_dict["subGT"]], params_dict["pInv"], params_dict["gamma"]
	new_orig_ll, tree_str = get_raxml_optimized_bl(tree_str_not_opt, dataset_path+MSA_PHYLIP_FILENAME, rates, pinv, alpha, freq)

	end1 = all_SPR(dataset_path, outpath, tree=tree_str)
	end2 = collect_features(dataset_path, relpath, outpath.format("prune"), outpath.format("rgft"))
	'''
	#end3 = os.system("python /groups/itay_mayrose/danaazouri/PhyAI/code/execute_programs/RF_learning.py -st exampleSphaero_{} -mt merged -sscore".format(relpath))
	os.rename(SUMMARY_FILES_DIR + "learning_all_moves_stepexampleSphaero_{}.csv".format(relpath), SUMMARY_FILES_DIR + "model_testing_exampleSphaero_st{}.csv".format(int(relpath.split("_")[-1][2:])))
	#learning_all_moves_stepexampleSphaero_best_pred_st1.csv
	# todo: run RF and save the best tree among the top 5 best predictions + its respective d_ll (to a file named relpath + '.txt')
	# todo: rerun this script with new tree + new orig ll
