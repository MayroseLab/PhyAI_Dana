import sys

sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")
from execute_programs.SPR_move import *
from execute_programs.RF_learning import *





def SPR_generator_forward(ds_path, outpath, tree, orig_ll):
	orig_msa_file = ds_path + MSA_PHYLIP_FILENAME
	suf = "bionj" if not RANDOM_TREE_DIRNAME in ds_path else "br"
	stats_filepath = ds_path + PHYML_STATS_FILENAME.format(suf) if ML_SOFTWARE_STARTING_TREE == 'phyml' else ds_path + RAXML_STATS_FILENAME

	t_orig = get_tree(ds_path, orig_msa_file, False) if not tree else PhyloTree(newick=tree, alignment=orig_msa_file, alg_format="iphylip", format=1)
	t_orig.get_tree_root().name = ROOTLIKE_NAME if not tree else ROOTLIKE_NAME + "_2"

	st = "1" if not tree else 'best_pred_st' + outpath.split("_")[-1][2:-4]
	OUTPUT_TREES_FILE = TREES_PER_DS.format(ds_path, st)
	with open(OUTPUT_TREES_FILE, "w", newline='') as fpw:
		csvwriter = csv.writer(fpw)
		csvwriter.writerow(['', 'prune_name', 'rgft_name', 'newick'])

	# first, copy msa file to memory and save it:
	msa_rampath = "/dev/shm/tmp" + ds_path.split(SEP)[-2] #  to be on the safe side (even though other processes shouldn't be able to access it)
	with open(orig_msa_file) as fpr:
		msa_str = fpr.read()

	try:
		with open(msa_rampath, "w") as fpw:
			fpw.write(msa_str)  # don't write the msa string to a variable (or write and release it)
		msa_str = ''

		params_dict = (parse_phyml_stats_output(None, stats_filepath)) if ML_SOFTWARE_STARTING_TREE == 'phyml' else parse_raxmlNG_output(stats_filepath)
		freq, rates, pinv, alpha = [params_dict["fA"], params_dict["fC"], params_dict["fG"], params_dict["fT"]], [params_dict["subAC"], params_dict["subAG"], params_dict["subAT"], params_dict["subCG"],params_dict["subCT"], params_dict["subGT"]], params_dict["pInv"], params_dict["gamma"]
		df = pd.DataFrame()
		for i, prune_node in enumerate(t_orig.iter_descendants("levelorder")):
			prune_name = prune_node.name
			nname, subtree1, subtree2 = prune_branch(t_orig, prune_name) # subtree1 is the pruned subtree. subtree2 is the remaining subtree
			with open(OUTPUT_TREES_FILE, "a", newline='') as fpa:
				csvwriter = csv.writer(fpa)
				csvwriter.writerow([str(i)+",0", prune_name, SUBTREE1, subtree1.write(format=1)])
				csvwriter.writerow([str(i)+",1", prune_name, SUBTREE2, subtree2.write(format=1)])

			for j, rgft_node in enumerate(subtree2.iter_descendants("levelorder")):
				ind = str(i) + "," + str(j)
				rgft_name = rgft_node.name
				if nname == rgft_name: # if the rgrft node is the one that was pruned
					continue
				rearr_tree_str = regraft_branch(subtree2, rgft_node, subtree1, rgft_name, nname).write(format=1)

				### save tree to file by using "append"
				with open(OUTPUT_TREES_FILE, "a", newline='') as fpa:
					csvwriter = csv.writer(fpa)
					csvwriter.writerow([ind, prune_name, rgft_name, rearr_tree_str])

				ll_rearr, rtime = call_raxml_mem(rearr_tree_str, msa_rampath, rates, pinv, alpha, freq)

				df.loc[ind, "prune_name"], df.loc[ind, "rgft_name"] = prune_name, rgft_name
				df.loc[ind, "prune_name"], df.loc[ind, "rgft_name"] = prune_name, rgft_name
				df.loc[ind, "time"] = rtime
				df.loc[ind, "ll"] = ll_rearr

		df["orig_ds_ll"] = orig_ll
		df.to_csv(outpath.format("prune"))
		df.to_csv(outpath.format("rgft"))

	except Exception as e:
		print('could not complete the all_SPR function on dataset:', dataset_path, '\nError message:')
		print(e)
		exit()
	finally:
		os.remove(msa_rampath)
	return



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
	dataset_name = dataset_path.split(SEP)[-2]
	analysis_st = int(relpath.split("_")[-1][2:]) + 1
	outpath = SUMMARY_PER_DS.format(dataset_path, "{}", "br", relpath)


	with open(dataset_path + relpath + '.txt', 'r') as fp:
		tree_str_not_opt = fp.read().strip()

	# update to string after bl optimization
	stats_filepath = dataset_path + PHYML_STATS_FILENAME.format("bionj")
	params_dict = (parse_phyml_stats_output(None,stats_filepath))
	freq, rates, pinv, alpha = [params_dict["fA"], params_dict["fC"], params_dict["fG"], params_dict["fT"]], [params_dict["subAC"], params_dict["subAG"], params_dict["subAT"], params_dict["subCG"], params_dict["subCT"],params_dict["subGT"]], params_dict["pInv"], params_dict["gamma"]
	new_orig_ll, tree_str = get_raxml_optimized_bl(tree_str_not_opt, dataset_path+MSA_PHYLIP_FILENAME, rates, pinv, alpha, freq)
	if not os.path.exists(SUMMARY_FILES_DIR + "model_testing_{}_st{}.csv".format(dataset_name, analysis_st)):
		end1 = SPR_generator_forward(dataset_path, outpath, tree_str, new_orig_ll)
		end2 = collect_features(dataset_path, relpath, outpath.format("prune"), outpath.format("rgft"))

		end3 = os.system("python {}execute_programs/RF_learning.py -st {}_{} -mt merged -sscore".format(CODE_PATH, dataset_name, relpath))
		end4 = os.rename(SUMMARY_FILES_DIR + "learning_all_moves_step{}_{}.csv".format(dataset_name, relpath), SUMMARY_FILES_DIR + "model_testing_{}_st{}.csv".format(dataset_name, analysis_st))
	if not os.path.exists(SUMMARY_FILES_DIR + DATA_WITH_PREDS.format('19_1_{}_st{}_ytransformed_exp'.format(dataset_name, analysis_st))):
		end5 = os.system('python {}execute_programs/RF_learning.py -st 1 -trans exp -mt merged -sscore -val {}_st{}'.format(CODE_PATH, dataset_name, analysis_st))

	df_with_preds = pd.read_csv(SUMMARY_FILES_DIR + DATA_WITH_PREDS.format('19_1_{}_st{}_ytransformed_exp'.format(dataset_name, analysis_st)))
	temp_df = df_with_preds.sort_values(by='pred', ascending=False).reset_index()
	dll, rank, prune, rgft = 0, None, None, None
	for i, row in temp_df.head(250).iterrows():
		true_dll = row[LABEL.format('prune')]
		if true_dll > dll:
			dll = true_dll
			prune = row["prune_name"]
			rgft = row["rgft_name"]
			rank = i + 1
			ll = row["ll"]

	if not prune:
		print("step{}:\nThe top predictions did not achieve likelihood improvement :(".format(analysis_st))
		print("log-lokelihood of the resulting final tree (best from previous step, bl-optimized): is: {}".format(new_orig_ll))
	else:
		print("step{}:\ndll folowing the #{} top prediction is: {}".format(analysis_st, rank, dll))
		print("(log-lokelihood: {})".format(ll))
		# locate this best tree in the respective dfr and save it to a file named relpath + '.txt'
		dfr = pd.read_csv(dataset_path + "newicks_step{}.csv".format(relpath))
		next_tree_str = dfr[(dfr["prune_name"] == prune) & (dfr["rgft_name"] == rgft)]["newick"].values[0]
		print(next_tree_str)
		with open (dataset_path + 'best_pred_st{}.txt'.format(analysis_st), 'w') as fp:
			fp.write(next_tree_str)

		# rerun this script for the next step
		#os.system("python {}execute_programs/SPR_forward.py -ds {} -trp best_pred_st{}".format(CODE_PATH, dataset_path, analysis_st))
