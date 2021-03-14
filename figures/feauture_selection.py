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
from matplotlib.pyplot import cm
from matplotlib import colors


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
	#plt.savefig(SUMMARY_FILES_DIR + "Fig2.tif", dpi=300)
	#plt.savefig(SUMMARY_FILES_DIR + "FigS1.tif", dpi=300)
	plt.show()
	
	return


def plot_main_results2_oneline(df):
	sns.set_context("paper", font_scale=1.2)
	#f, (axA1, axA2, axA3), (axB1, axB2, axB3) = plt.subplots(2,3, figsize=(15,12))
	# axes = [axA1, axA2, axA3, axB1, axB2, axB3]
	f, axes = plt.subplots(2,3, figsize=(15,10))
	(axA1, axA2, axA3) = axes[0][0], axes[0][1], axes[0][2]
	(axB1, axB2, axB3) = axes[1][0], axes[1][1], axes[1][2]
	my_colors = []
	for i, ax in enumerate([axA1, axA2, axA3]):
		my_colors.append(next(palette))
		sns.distplot(df[SCORES_LST[i]].values, ax=ax, color=my_colors[i], kde=False, label="big",bins=25)  # , norm_hist=True)

	ypercA1 = [int(100 * y / len(df)) for y in axA1.get_yticks()]
	#yperc2 = [int(100 * y / len(df)) for y in ax2.get_yticks()]
	#yperc3 = [int(100 * y / len(df)) for y in ax3.get_yticks()]
	
	axA1.set(ylabel="% Empirical datasets", yticklabels=ypercA1,xlim=(0, 1))
	axA2.set(xlim=(0, 100))
	axA3.set(xlim=(0, 100))
	axA2.set_yticklabels(axA2.get_yticklabels(), size=30)
	axA3.set_yticklabels(axA3.get_yticklabels(), size=30)
	
	axA1.text(0, 1.08, "a$_{1}$", transform=axA1.transAxes, fontsize=18, fontweight='bold', va='top', ha='right')
	axA2.text(0, 1.08, "a$_{2}$", transform=axA2.transAxes, fontsize=18, fontweight='bold', va='top', ha='right')
	axA3.text(0, 1.08, "a$_{3}$", transform=axA3.transAxes, fontsize=18, fontweight='bold', va='top', ha='right')
	
	###################
	df = pd.read_csv("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/summary_files/results_saturation/scores_per_ds_20_1_validation_set_4200_ytransformed_exp10_cp.csv")
	for i, ax in enumerate([axB1, axB2, axB3]):
		sns.distplot(df[SCORES_LST[i]].values, ax=ax, color=my_colors[i], kde=False, label="big",bins=25)  # , norm_hist=True)
	ypercB1 = [int(100 * y / len(df)) for y in axB1.get_yticks()]
	
	axB1.set(xlabel='Spearman correlation ({})'.format(r'$\rho$'), ylabel="% Empirical datasets", yticklabels=ypercB1, xlim=(0, 1))
	axB2.set(xlabel="Empirically best ranking percentile", xlim=(0, 100))
	axB3.set(xlabel="Best predicted ranking percentile", xlim=(0, 100))
	axB2.set_yticklabels(axB2.get_yticklabels(), size=30)
	axB3.set_yticklabels(axB3.get_yticklabels(), size=30)
	
	axB1.text(0, 1.08, "b$_{1}$", transform=axB1.transAxes, fontsize=18, fontweight='bold', va='top', ha='right')
	axB2.text(0, 1.08, "b$_{2}$", transform=axB2.transAxes, fontsize=18, fontweight='bold', va='top', ha='right')
	axB3.text(0, 1.08, "b$_{3}$", transform=axB3.transAxes, fontsize=18, fontweight='bold', va='top', ha='right')
	###################
	
	f.tight_layout()
	#plt.savefig(SUMMARY_FILES_DIR + "FigS4.tif", dpi=300)
	#plt.savefig(SUMMARY_FILES_DIR + "Fig2.tif", dpi=300)
	plt.show()

	return



removed_featured_dict = {19: 'Sum of branches between pruning and regrafting',
		   18: 'Number of branches between pruning and regrafting',
		   16: 'The length of the pruned branch',
		   15: 'The longest branch of the subtree in Fig. 1b',
		   17: 'The sum of branches in the starting tree',
		   14: 'The sum of branches of the subtree in Fig. 1b',
		   12: 'The sum of branches of the subtree in Fig. 1c',
		   9: 'The length of the regrafted branch',
		   13: 'The longest branch in the starting tree',
		   11: 'The length of the branch formed due to pruning',
		   8: 'The longest branch of the subtree in Fig. 1c',
		   7: 'The sum of branches of the subtree in Fig. 1c2',
		   4: 'The longest branch of the subtree in Fig. 1c2',
		   10: 'The longest branch of the subtree in Fig. 1c1',
		   5: 'The sum of branches of the subtree in Fig. 1c1',
		   6: 'The number of leaves in the subtree in Fig. 1c',
		   3: 'The number of leaves in the subtree in Fig. 1b',
		   2: 'The number of leaves in the subtree in Fig. 1c2',
		   1: 'The number of leaves in the subtree Fig. 1c1'}


def scores_feature_selection(df):
	fig = plt.figure(figsize=(14, 10))
	colors_list = cm.YlGnBu(np.linspace(0, 0.8, 19))[::-1]
	sns.boxplot(x=N_FEATURES_COL, y=SCORES_LST[0], data=df, order=[19-i for i in range(19)], showfliers=False, saturation=0.6, linewidth=0.6, palette=colors_list)
	plt.xlabel('')
	plt.ylabel('Spearman correlation ({})'.format(r'$\rho$'), size=18)

	rows_text_nested_lst, colors_nested_lst = [], []
	for ind in range(19,0,-1):
		rows_text_nested_lst.append(['+' if ind>=x else '' for x in range(1,20)])
		colors_nested_lst.append([colors_list[x-1] if ind>=x else 'black' for x in range(1,20)])

	the_table = plt.table(cellText=rows_text_nested_lst, rowLoc='left',
						  rowLabels=[removed_featured_dict[19-i] for i in range(19)],
						  colLabels=[str(19-i) for i in range(19)],
						  colLoc='center',cellLoc = 'center', cellColours=colors_nested_lst, # rowColours=['whitesmoke' for x in range(1,21)],
						  loc='bottom')

	the_table.auto_set_font_size(False)
	the_table.set_fontsize(14)
	the_table.scale(1., 1.7)

	plt.xticks([])
	fig.subplots_adjust(left=0.2 ,bottom=0.32)
	plt.text(7, -0.47, 'Features in the model', size=18)

	plt.tight_layout()
	plt.savefig(dirpath + "feature_selection_larger_mat.tif", dpi=300)
	#plt.show()






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
	'''
	##dirpath += 'results_feature_selection/'
	#df = concat_n_features(dirpath, max_n_features=19)
	#df.to_csv(dirpath + 'temp.csv')
	scores_feature_selection(pd.read_csv(dirpath + 'temp.csv'))
	#exit()
	df = pd.read_csv(dirpath + 'temp.csv')
	a8 = df.loc[df[N_FEATURES_COL] == 8, SCORES_LST[0]].dropna().values
	a9 = df.loc[df[N_FEATURES_COL] == 9, SCORES_LST[0]].dropna().values
	a10 = df.loc[df[N_FEATURES_COL] == 10, SCORES_LST[0]].dropna().values
	a11 = df.loc[df[N_FEATURES_COL] == 11, SCORES_LST[0]].dropna().values
	a12 = df.loc[df[N_FEATURES_COL] == 12, SCORES_LST[0]].dropna().values
	a13 = df.loc[df[N_FEATURES_COL] == 13, SCORES_LST[0]].dropna().values
	a14 = df.loc[df[N_FEATURES_COL] == 14, SCORES_LST[0]].dropna().values
	a15 = df.loc[df[N_FEATURES_COL] == 15, SCORES_LST[0]].dropna().values
	a16 = df.loc[df[N_FEATURES_COL] == 16, SCORES_LST[0]].dropna().values
	a17 = df.loc[df[N_FEATURES_COL] == 17, SCORES_LST[0]].dropna().values
	a18 = df.loc[df[N_FEATURES_COL] == 18, SCORES_LST[0]].dropna().values
	a19 = df.loc[df[N_FEATURES_COL] == 19, SCORES_LST[0]].dropna().values
	import scipy.stats as stats

	F, p = stats.f_oneway(a19,a18,a17,a16, a15)
	print(p)
	F, p = stats.f_oneway(a19, a14, a13, a12, a11, a10, a9, a8)
	print(p)
	F, p = stats.f_oneway(a19, a14)
	print(p)
	F, p = stats.f_oneway(a19, a13)
	print(p)
	F, p = stats.f_oneway(a19, a12)
	print(p)
	F, p = stats.f_oneway(a19, a11)
	print(p)
	F, p = stats.f_oneway(a19, a10)
	print(p)
	F, p = stats.f_oneway(a19, a9)
	print(p)
	F, p = stats.f_oneway(a19, a8)
	print(p)
	exit()
	#'''

	df = pd.read_csv(dirpath + 'v4_20features/' + SCORES_PER_DS.format("20_1_4200_ytransformed_exp"))
	#df = pd.read_csv("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/summary_files/results_saturation/scores_per_ds_20_1_validation_set_4200_ytransformed_exp10_cp.csv")
	#df = pd.read_csv("/groups/itay_mayrose/danaazouri/PhyAI/submission_data/summary_files/scores_per_ds_19_1_validation_set_ytransformed_exp.csv")  # JC
	plot_main_results2_oneline(df)
	#plot_main_results2(df)

	
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