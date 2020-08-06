import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")
import warnings
warnings.filterwarnings("ignore")			# TEMP

from defs import *

import seaborn as sns
sns.set_style("white")
N_BINS_y = 3
N_BINS = 3

SAMPLE_SIZE_METRICS = ["ntaxa"] #,nsamples]
CONTENT_TYPE = "best predicted in true" # best predicted in true | ntaxa_needed
CONTENT_SETTINGS = {"best predicted in true": {"ylim": [[0,30], [0,30], [0,30]], "ylabel": "best pridicted rank"}}
					#,"nsamples_needed": {"ylim": [[-0.01,20], [-0.01,20], [-0.01,20]], "ylabel": "n samples needed"},}
dirpath = SUMMARY_FILES_DIR if platform.system() == 'Linux' else DATA_PATH



def accXsize_boxplot(res_vec):
	plt.gcf().subplots_adjust(bottom=0.15, left=0.30)
	f, axarr = plt.subplots(N_BINS_y, 1, figsize=(6, 6))
	fontsize = 9

	checked_criteria = ["rank score"]
	data_summary = pd.read_csv(dirpath + CHOSEN_DATASETS_FILENAME)
	print(len(data_summary))
	print(len(res_vec))
	data_summary[checked_criteria[0]] = res_vec

	# binning
	binned_category = data_summary["nchars"]
	data_summary['nchars_cat'], bins = pd.qcut(binned_category, q=N_BINS, retbins=True, precision=2)
	bins = [int(bin) for bin in bins]
	bins[0] -= 1

	binned_category_y = data_summary["ntaxa"]
	data_summary['ntaxa_cat'], bins_y = pd.qcut(binned_category_y, q=N_BINS_y, retbins=True, precision=2, duplicates='drop')
	bins_y = [int(bin) for bin in bins_y]
	#bins_y[0] -= 1

	for j in range(N_BINS_y):
		#filtered_idx_y = (bins_y[j] < binned_category_y) & (binned_category_y <= bins_y[j + 1])
		filtered_idx_y = (bins_y[j] == binned_category_y)
		current_ntaxa_df = copy.deepcopy(data_summary)
		current_ntaxa_df = current_ntaxa_df.ix[filtered_idx_y]

		for criterion in checked_criteria:

			for i in range(N_BINS):
				current_ntaxa_df.ix[(bins[i] < binned_category) &
									(binned_category <= bins[i + 1]), "nchars_bin"] = \
					"({0}, {1}]".format(bins[i], bins[i + 1])

		whole_df = pd.melt(current_ntaxa_df, id_vars=['nchars_bin'], value_vars=checked_criteria,
						   var_name="criterion", value_name=CONTENT_SETTINGS[CONTENT_TYPE]["ylabel"])

		bp = sns.violinplot(x="nchars_bin", y=CONTENT_SETTINGS[CONTENT_TYPE]["ylabel"], inner="point", scale="width",
						 data=whole_df, hue="criterion", color="dodgerblue", width=0.8, ax=axarr[j], # , palette="Set3"
						 sym='k.', showmeans=False, meanline=True, showfliers=True, linewidth=0.5,
						 order=["({0}, {1}]".format(bins[i], bins[i + 1]) for i in range(N_BINS)],
						 )

		axarr[j].set_ylim(CONTENT_SETTINGS[CONTENT_TYPE]["ylim"][j])
		bins[0] += 1

		axarr[j].legend_.remove()

		if j != N_BINS - 1:
			axarr[j].set_xlabel("")
			axarr[j].set_xticks([])
		else:
			axarr[j].set_xticklabels(axarr[j].get_xticklabels(), fontsize=fontsize)

		#axarr[j].set_ylabel("{}-{}".format(bins_y[j] + 1, bins_y[j + 1]), rotation=270, labelpad=10)
		axarr[j].set_ylabel(bins_y[j], rotation=270, labelpad=10)
		axarr[j].yaxis.set_label_position("right")

	axarr[1].text(-0.65, 22, CONTENT_SETTINGS[CONTENT_TYPE]["ylabel"], rotation=90, fontsize=fontsize + 2)
	axarr[1].text(2.65, 17, "# Taxa", rotation=270, fontsize=fontsize + 2)

	sns.despine(offset=0, left=False)


	axarr[0].legend(loc="best", bbox_to_anchor=(1.18, 1.0), fontsize=fontsize, frameon=False) #, bbox_to_anchor=(1.8, 1.0)
	plt.xlabel("MSA length", fontsize=fontsize + 2, labelpad=20)

	plt.tight_layout(h_pad=0.05)
	f.set_size_inches(12, 9, forward=True)

	plt.savefig("samplesize.png", dpi=300, bbox_inches="tight")
	plt.show()




if __name__ == '__main__':
	# test. not real
	accXsize_boxplot(np.concatenate((np.arange(970), np.asarray([2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 5, 6, 7, 8, 16, 25, 27, 27, 53]))), 7)