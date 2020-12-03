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

		p = Popen([RAXML_NG_SCRIPT, '--evaluate', '--msa', msa_file,'--threads', '2', '--opt-branches', 'on', '--opt-model', 'off', '--model', model_line_params, '--nofiles', '--tree', tree_rampath], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
		raxml_stdout = p.communicate()[0]
		raxml_output = raxml_stdout.decode()
		print("\n"+raxml_output+"\n")
		exit()

		# todo: parse to return also tree str
		res_dict = parse_raxmlNG_content(raxml_output)
		ll = res_dict['ll']
		rtime = res_dict['time']

	except Exception as e:
		print(msa_file.split(SEP)[-1][3:])
		print(e)
		exit()
	finally:
		os.remove(tree_rampath)

	# todo: return also (optimized) tree str
	return ll, rtime






if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='perform all SPR moves')
	parser.add_argument('--dataset_path', '-ds', default=None)
	parser.add_argument('--tree_relpath', '-trp', required=True)	# e.g., best_pred_st13
	args = parser.parse_args()

	dataset_path = args.dataset_path
	outpath = SUMMARY_PER_DS.format(dataset_path, "{}", "br", args.tree_relpath)

	res = all_SPR(dataset_path, outpath, tree=None, rewrite_phylip=args.rewrite_in_phylip)

	with open(dataset_path + args.tree_relpath + '.txt', 'r') as fp:
		tree_str_not_opt = fp.read().strip()

	# update to string after bl optimization
	stats_filepath = dataset_path + PHYML_STATS_FILENAME.format("bionj")
	params_dict = (parse_phyml_stats_output(None,stats_filepath))
	freq, rates, pinv, alpha = [params_dict["fA"], params_dict["fC"], params_dict["fG"], params_dict["fT"]], [params_dict["subAC"], params_dict["subAG"], params_dict["subAT"], params_dict["subCG"], params_dict["subCT"],params_dict["subGT"]], params_dict["pInv"], params_dict["gamma"]
	tree_str = get_raxml_optimized_bl(tree_str_not_opt, dataset_path+MSA_PHYLIP_FILENAME, rates, pinv, alpha, freq)
	res = all_SPR(dataset_path, outpath, tree=tree_str, rewrite_phylip=args.rewrite_in_phylip)

	collect_features(dataset_path, args.tree_relpath, outpath.format("prune"), outpath.format("rgft"))
	# todo: run RF and extract the best 5 best predictions
	# todo: save the ~best predicted tree + additional inforamtion
	# todo: rerun this script with new tree + new orig ll
