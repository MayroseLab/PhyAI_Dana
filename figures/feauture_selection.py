import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")
import warnings
warnings.filterwarnings("ignore")			# TEMP
import itertools
import seaborn as sns

from defs import *
N_FEATURES_COL = "n_features"

sns.set_style("white")
palette = itertools.cycle(sns.color_palette('colorblind'))

import matplotlib.gridspec as gridspec


def plot_main_results(df):
	f = plt.figure()
	gs0 = gridspec.GridSpec(1, 3)
	gs00 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[0])
	gs01 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[1])
	gs02 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[2])

	ax1, ax2, ax3, ax4 = f.add_subplot(gs00[:]), f.add_subplot(gs01[:,:1]), f.add_subplot(gs01[:,1:2]),f.add_subplot(gs02[:])

	ylims = [(0, 1), (1, 16000), (1, 16000), (0, 100)]
	axes = [ax1,ax2,ax3,ax4]
	for i,ax in enumerate(axes):
		sns.violinplot(y=SCORES_LST[i], data=df, ax=ax, color=next(palette))
		ax.set_ylim(ylims[i])


	ax1.set(xlabel="Spearman correlation", ylabel="Score (0 to 1)")
	ax2.set(xlabel="Best predicted\nranking", ylabel="Score (1 to 100)")
	ax3.set(xlabel="Empirically best\nranking", ylabel="")
	ax4.set(xlabel="Required evaluation set", ylabel="Score (0% to 100%)")

	#for i,ax in enumerate(axes):
	#	ax.set_xticklabels(ax.get_xticklabels(), size=30)
	#	ax.set_yticklabels(ax.get_yticklabels(), size=30)

	f.tight_layout()
	plt.show()
	#plt.savefig(dirpath + "results.png")


def plot_main_results2(df):
	f = plt.figure()
	sns.set_context("paper", font_scale=1.2)
	gs0 = gridspec.GridSpec(2, 2)
	gs00 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[0,:])
	gs01 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[1,:])

	ax1, ax2, ax3 = f.add_subplot(gs00[:]), f.add_subplot(gs01[:, :1]), f.add_subplot(gs01[:, 1:2])

	#ylims = [(0, 1), (1, 16000), (1, 16000)]
	axes = [ax1, ax2, ax3]
	for i, ax in enumerate(axes):
		sns.distplot(df[SCORES_LST[i]].values, ax=ax, color=next(palette), kde=False, label="big", bins=25)#, norm_hist=True)
	
	yperc1 = [int(100*y / len(df)) for y in ax1.get_yticks()]
	yperc2 = [int(100 * y / len(df)) for y in ax2.get_yticks()]

	ax1.set(xlabel='Spearman correlation ({})'.format(r'$\rho$'), ylabel="% Empirical datasets", yticklabels=yperc1, xlim=(0,1))
	ax2.set(xlabel="Empirically best ranking percentile", ylabel="% Empirical datasets", yticklabels=yperc2, xlim=(0,100))
	ax3.set(xlabel="Best predicted ranking percentile", xlim=(0,100))
	ax3.set_yticklabels(ax3.get_yticklabels(), size=30)
	
	ax1.text(0, 1.15, "a", transform=ax1.transAxes, fontsize=18, fontweight='bold', va='top', ha='right')
	ax2.text(0, 1.15, "b", transform=ax2.transAxes, fontsize=18, fontweight='bold', va='top', ha='right')
	ax3.text(0, 1.15, "c", transform=ax3.transAxes, fontsize=18, fontweight='bold', va='top', ha='right')

	f.tight_layout()
	f.set_size_inches(7, 7, forward=True)
	plt.savefig("C:\\Users\\ItayMNB3\\Dropbox\\PhyloAI\\PhyAI_writing\\to_submit\\" + "Fig2.tif", dpi=300)
	#plt.savefig("C:\\Users\\ItayMNB3\\Dropbox\\PhyloAI\\PhyAI_writing\\figures\\" + "FigS1.tif", dpi=300)
	#plt.show()
	
	return



def scores_feature_selection(df):
	'''
	fig, axarr = plt.subplots(2, 2)
	axarrs_locs = [(0, 0), (0, 1), (1, 0)]
	ylims = [(0.0, 1), (1, 100), (1, 100)]
	for i,loc in enumerate(axarrs_locs):
		color = next(palette)
		ax = axarr[axarrs_locs[i]]
		ax.set_ylim(ylims[i])
		sns.boxplot(x=N_FEATURES_COL, y=SCORES_LST[i], data=df, ax=ax, color=color, showfliers=False)
		ax.set_xticklabels(ax.get_xticklabels(), rotation=30)
	'''
	removed_featured_dict = {20: 'The full set of features',
			   16: 'The number of branches in the path between the regrafting and the pruning branches',
			   19: 'The estimated total branch lengths of the resulting tree',
			   17: 'The length of the branch that was being pruned',
			   18: 'The length of the longest branch in the subtree in Fig. 1b',
			   15: 'The sum of branches in the starting tree',
			   12: 'The sum of branches in the subtree in Fig. 1b',
			   10: 'The sum of branches in the subtree in Fig. 1c',
			   8: 'The length of the branch that was being regrafted',
			   13: 'The length of the longest branch in the starting tree',
			   14: 'The approximated length of the newly formed branch formed due to regrafting',
			   9: 'The length of the longest branch in the subtree in Fig. 1c',
			   6: 'The sum of branches in the subtree in Fig. 1c2',
			   5: 'The length of the longest branch in the subtree in Fig. 1c2',
			   11: 'The length of the longest branch in the subtree in Fig. 1c1',
			   4: 'The sum of branches in the subtree in Fig. 1c1',
			   7: 'The number of leaves in the subtree in Fig. 1c',
			   3: 'The number of leaves in the subtree in Fig. 1b',
			   2: 'The number of leaves in the subtree in Fig. 1c2',
		       1: 'The number of leaves in the subtree Fig. 1c1'}



	fig, axarr = plt.subplots()
	sns.boxplot(x=N_FEATURES_COL, y=SCORES_LST[0], data=df, showfliers=False, saturation=0.6, linewidth=0.6)
	#plt.xlabel('Number of features in the model')
	plt.xlabel('')
	plt.ylabel('Spearman correlation ({})'.format(r'$\rho$'))

	#axarr.set_xticklabels([removed_featured_dict[i] + " ({})".format(i) for i in range(1,21)], rotation=90)
	# Add a table at the bottom of the axes
	rows_text_nested_lst = []
	for ind in range(1,21):
		rows_text_nested_lst.append([1 for x in range(ind)] )

	'''
	the_table = plt.table(cellText=['+' for i in range(1,21)],
						  rowLabels=[removed_featured_dict[i] for i in range(1,21)],
						  #rowColours=colors,
						  colLabels=[str(i) for i in range(1,21)],
						  loc='bottom')

	# Adjust layout to make room for the table:
	#plt.subplots_adjust(left=0.2, bottom=0.2)
	'''

	import matplotlib.cm as cm


	fig = plt.figure()
	ax1 = fig.add_subplot(121)
	# Bilinear interpolation - this will look blurry
	ax1.imshow(rows_text_nested_lst, interpolation='bilinear', cmap=cm.Greys_r)

	#ax2 = fig.add_subplot(122)
	# 'nearest' interpolation - faithful but blocky
	#ax2.imshow(rows_text_nested_lst, interpolation='nearest', cmap=cm.Greys_r)

	plt.tight_layout()
	#fig.set_size_inches(7, 4, forward=True)
	#plt.savefig(dirpath + "feature_selection_with_names.tif", dpi=300)
	plt.show()


def concat_n_features(dirpath, max_n_features):
	'''
	pre-merge iterations results to one csv with "n_features" col
	'''
	df = pd.read_csv(dirpath + SCORES_PER_DS.format(str(max_n_features) + "_1_ytransformed_exp"))
	df[N_FEATURES_COL] = max_n_features
	for i in range(max_n_features-1,0,-1):
		csv = dirpath + SCORES_PER_DS.format(str(i) + "_1_ytransformed_exp")
		if not os.path.exists(csv):
			break
		df2 = pd.read_csv(csv)
		df2[N_FEATURES_COL] = i
		df = pd.concat([df,df2])

	return df


	


if __name__ == '__main__':
	dirpath = SUMMARY_FILES_DIR if platform.system() == 'Linux' else DATA_PATH
	dirpath += 'results_feature_selection/'
	df = concat_n_features(dirpath, max_n_features=20)
	scores_feature_selection(df)
	
	#df = pd.read_csv(dirpath + SCORES_PER_DS.format("28_1_ytransformed_exp"))
	#df = pd.read_csv(dirpath + SCORES_PER_DS.format("26_2_ytransformed_exp"))
	#df = pd.read_csv(dirpath + SCORES_PER_DS.format("26_1st_on_2nd_ytransformed_exp"))

	#df = pd.read_csv(dirpath + SCORES_PER_DS.format("28_1_ytransformed_exp_minus"))
	##df = pd.read_csv(dirpath + SCORES_PER_DS.format("20_1_3700_ytransformed_exp"))
	#df = pd.read_csv(dirpath + SCORES_PER_DS.format("26_1_validation_set_ytransformed_exp"))
	##plot_main_results2(df)
	
	
	'''
	temp_dir = "C:\\Users\\ItayMNB3\\Desktop\\"
	# temp_dir = "D:\\Users\\Administrator\\Desktop\\"
	df = pd.read_csv(temp_dir + "with_preds_merged_26_1_validation_set_ytransformed_exp.csv")
	cnt1, cnt_groups = 0, 0
	
	grouped_df_by_ds = df.groupby(FEATURES["group_id"], sort=False)
	for group_id, df_by_ds in grouped_df_by_ds:
		cnt_groups +=1
		best_pred_ix = df_by_ds["pred"].idxmax()  # changed min to max!

		reverse_transformed = (df_by_ds["orig_ds_ll"].values[0]) * (
			np.log2(df_by_ds.loc[df_by_ds["Unnamed: 0"] == best_pred_ix]["pred"] - 1)).values
		#print(reverse_transformed)
		#print(df_by_ds["orig_ds_ll"].values[0])
		#print(df_by_ds["d_ll_prune"].values[0])
		#exit()
		if reverse_transformed > 0:
			cnt1 += 1

	print(cnt1)
	print(100*cnt1/cnt_groups)
	print(100*len(df[df["d_ll_prune"] > 0])/len(df))
	'''