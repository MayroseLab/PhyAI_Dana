import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")
import warnings
warnings.filterwarnings("ignore")			# TEMP
import itertools
import seaborn as sns

from defs import *
from utils.tree_functions import *

sns.set_style("white")
#fig, axarr = plt.subplots(2, 2)
palette = itertools.cycle(sns.color_palette('colorblind'))
from scipy.stats import pearsonr
from Bio import AlignIO


def count_gaps_proportion(msa):
	alignment = AlignIO.read(msa, PHYLIP_FORMAT)
	avg_gappiness = []
	for record in alignment:
		avg_gappiness.append(record.seq.count("-"))
	avg_gappiness_per_taxa = np.mean(avg_gappiness)
	gappiness = round(100*avg_gappiness_per_taxa / alignment.get_alignment_length(), 2)

	return gappiness


def get_node_height(tree_path):
	root_node = Tree(tree_path, format=1).get_tree_root()
	heights = []
	for i, leaf in enumerate(root_node.iter_leaves()):
		heights.append(root_node.get_distance(leaf))
	var = round(np.var(heights), 5)

	return var


def calc_empirical_features():
	df = pd.read_csv('/groups/itay_mayrose/danaazouri/PhyAI/validation_set2/summary_files/sampled_datasets.csv').reset_index(drop=True)
	df2 = pd.read_csv("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/summary_files/results_saturation/scores_per_ds_20_1_validation_set_4200_ytransformed_exp_orig.csv").reset_index(drop=True)
	df2 = df2[df2["path"] != "/groups/itay_mayrose/danaazouri/PhyAI/DBset2/data/training_datasets/41157/"]
	for index, row in df2.iterrows():
		path = row["path"]
		sum_ds_df = SUMMARY_PER_DS.format(path, "prune", "br", "1")
		if os.path.exists(sum_ds_df):
			df2.loc[index, "ntaxa"] = df.loc[df["path"] == path, "ntaxa"].values[0]
			df2.loc[index, "nchars"] = df.loc[df["path"] == path, "nchars"].values[0]
			df2.loc[index, "tbl"] = get_total_branch_lengths(path + PHYML_TREE_FILENAME.format('bionj'))
			df2.loc[index, "ll"] = pd.read_csv(sum_ds_df).loc[0, "orig_ds_ll"]
			df2.loc[index, "gaps"] = count_gaps_proportion(row["path"] + MSA_PHYLIP_FILENAME)
			df2.loc[index, "theight_var"] = get_node_height(row["path"] + PHYML_TREE_FILENAME.format('bionj'))

	df2.to_csv(SUMMARY_FILES_DIR + SCORES_PER_DS.format("validation_with_atts"))



def plot_distributions(df):
	axarrs_locs = [(0, 0), (0, 1), (1, 0), (1, 1)]
	#xlims = [(0, 75), (0, 6000), (0, 50), (-50000, 0)]
	cols_lst = ["ntaxa", "nchars", "tbl", "ll"]
	for i, loc in enumerate(axarrs_locs):
		color = next(palette)
		ax = axarr[axarrs_locs[i]]
		#ax = axarr[i]
		#ax.set_xlim(xlims[i])
		'''
		for df_id in ["2_", "val_"]:
			df = pd.read_csv(DATA_PATH + df_id + CHOSEN_DATASETS_FILENAME, usecols=cols_lst)
			df_name = "TreeBASE set" if df_id == "2_" else "validation set"
			linestyle = "solid" if df_id == "2_" else 'dashed'
			df[cols_lst[i]].plot.kde(color=color, linestyle=linestyle, title=cols_lst[i], label=df_name, ax=ax)
			ax.legend(loc='best', frameon=True)
		'''
		q=5
		df[cols_lst[i]] = pd.qcut(df[cols_lst[i]], q=q) #, labels=np.arange(1,q+1))  #+'_attribute_category'
		df = df.round(2)
		#df2 = pd.read_csv(DATA_PATH + "2_" + CHOSEN_DATASETS_FILENAME)
		#df = pd.concat([df,df2])
		sns.boxplot(x=cols_lst[i], y=SCORES_LST[0], data=df,ax=ax, showfliers=False, color=color) #, hue='Database')  # ,color=color
		ax.legend(loc='best', frameon=True)
		
	fig.tight_layout()
	#plt.savefig(DIRPATH + 'dbs_attributes.png')
	plt.show()
	
	'''
	for df_id in ["2_", "val_"]:
		df = pd.read_csv(DATA_PATH + df_id + CHOSEN_DATASETS_FILENAME, usecols=["ntaxa","nchars","tbl","ll"])
		df_name = "TreeBASE set" if df_id=="2_" else "validation set"
		linestyle = "solid" if df_id=="2_" else 'dashed'
		ax = df.plot.kde(subplots=True,linestyle=linestyle,title=["ntaxa","nchars","tbl","ll"], label=False, layout = (2,2), ax=ax)   #label=df_name,

		# df.hist(color=color,linestyle=linestyle) #, ax=ax)
		'''


import matplotlib.gridspec as gridspec
from scipy.stats import *
import matplotlib.patches as mpatches
#import SeabornFig2Grid as sfg


def corr_plot(df):
	sns.set_context("paper", font_scale=1.5)
	palette = itertools.cycle(sns.color_palette('colorblind'))

	ax1 = sns.jointplot(x="ntaxa", y=SCORES_LST[0], data=df, kind='reg', stat_func=pearsonr, line_kws={'color':'black'}, color=next(palette),xlim=(5, 70.3), ylim=(0,1))#, ax=ax1)
	plt.xlabel('Number of taxa')
	plt.ylabel('Spearman correlation ({})'.format(r'$\rho$'))
	#stats_patch = mpatches.Patch(color='white', label='$r^2$ = 0.0002;  $pval$ = 0.34')#, contains=False)
	stats_patch = mpatches.Patch(color='white', label='$r^2$ = 0.009;  $pval$ = $7.8x10^-$$^1$$^0$')#, contains=False)
	plt.legend(handles=[stats_patch])
	plt.text(1.4, 1.2, "a", fontsize=20, fontweight='bold', va='top', ha='right')
	plt.tight_layout()
	plt.show()

	'''
	df = df[df["tbl"] < 30]
	ax2 = sns.jointplot(x="tbl", y=SCORES_LST[0], data=df, kind='reg', stat_func=pearsonr, line_kws={'color':'black'}, color=next(palette),xlim=(-0.3, 31), ylim=(0,1))#, ax=ax2)
	plt.xlabel('Total branch lengths')
	plt.ylabel("")
	# plt.xscale("log")
	#stats_patch = mpatches.Patch(color='white', label='$r^2$ = 0.007;  $pval$ = $1.1x10^-$$^8$')#, contains=False)
	stats_patch = mpatches.Patch(color='white', label='$r^2$ = 0.025;  $pval$ = $3.1x10^-$$^2$$^4$')  # , contains=False)
	plt.legend(handles=[stats_patch])
	plt.text(-1.9, 1.2, "b", fontsize=20, fontweight='bold', va='top', ha='right')
	plt.tight_layout()
	plt.show()
	'''
	f, axes = plt.subplots(3, 1, sharex=True, sharey=True)
	palette = itertools.cycle(sns.color_palette('colorblind'))
	df = df[df["tbl"] < 30]
	df1, df2, df3 = df[df["ntaxa"] <= 27], df[(df["ntaxa"] > 27) & (df["ntaxa"] <= 48)], df[df["ntaxa"] > 48]
	for i, df_i in enumerate([df1, df2, df3]):
		sns.regplot(x="tbl", y=SCORES_LST[0], data=df_i, line_kws={'color': 'black'}, color=next(palette), ax=axes[i])
		axes[i].yaxis.set_label_position("right")
		ilow, ihigh = (7, 27) if i == 0 else ("# Taxa \n\n28", 48) if i == 1 else (49, 70)
		axes[i].set_ylabel("{}-{}".format(ilow, ihigh), rotation=-90, labelpad=18) if i !=1 else axes[i].set_ylabel("{}-{}".format(ilow, ihigh), rotation=-90, labelpad=50)
		axes[i].set_xlabel('Total branch lengths') if i == 2 else axes[i].set_xlabel('')
		df_i = df_i[['tbl', SCORES_LST[0]]].dropna()
		r,p = pearsonr(df_i['tbl'].values, df_i[SCORES_LST[0]].values)
		stats_patch = mpatches.Patch(color='white',label='$r^2$ = {};  $pval$ = ${}x10^-$$^{}$$^{}$'.format((str(np.square(r)))[:5], str(p)[:4], str(p)[-2], str(p)[-1]))
		axes[i].legend(handles=[stats_patch], loc='lower right')
	axes[0].text(-1.4, 1.5, "b", fontsize=20, fontweight='bold', va='top', ha='right')
	plt.xlim(-0.3, 31)
	plt.ylim(0, 1)
	plt.tight_layout()
	plt.show()

	df_val = pd.read_csv(dirpath + SCORES_PER_DS.format("20_1_validation_set_ytransformed_exp"))
	df_val["Database"] = df_val["Database"].replace(["orthoMam"], "OrthoMaM")
	df_val["Set"] = "validation"
	df["Set"] = "training"
	df_inc_val = pd.concat([df, df_val])
	ax3 = sns.boxplot(x="Database", y=SCORES_LST[0], data=df_inc_val, showfliers=False,order=["ploiDB","protDBs", "selectome", "TreeBASE", "OrthoMaM", "PANDIT"], color="lightgrey")#, ax=ax3)
	plt.ylim(0,1)
	plt.ylabel("")
	plt.setp(ax3.get_xticklabels(), rotation=20)
	plt.text(-0.7, 1.2, "c", fontsize=20, fontweight='bold', va='top', ha='right')
	plt.tight_layout()
	plt.show()


	''''######################
	gpd_df = df_inc_val.groupby('Database')
	for i, df_db in gpd_df:
		print(df_db['Database'].values[0])
		print(df_db[SCORES_LST[0]].values.mean())
	######################'''




def corr_plot_more_atts(df):
	sns.set_context("paper", font_scale=1.5)
	palette = itertools.cycle(sns.color_palette('colorblind'))
	next(palette)

	ax1 = sns.jointplot(x="nchars", y=SCORES_LST[0], data=df, kind='reg', stat_func=pearsonr, line_kws={'color':'black'}, color=next(palette), ylim=(0,1))#, ax=ax1)
	plt.xlabel('Alignment length')
	plt.ylabel('Spearman correlation ({})'.format(r'$\rho$'))
	stats_patch = mpatches.Patch(color='white', label='$r^2$ = 0.01;  $pval$ = $5.2x10^-$$^1$$^1$')#, contains=False)
	plt.legend(handles=[stats_patch])
	plt.text(-3.4, 1.2, "d", fontsize=20, fontweight='bold', va='top', ha='right')
	plt.tight_layout()
	plt.show()

	ax1 = sns.jointplot(x="gaps", y=SCORES_LST[0], data=df, kind='reg', stat_func=pearsonr, line_kws={'color':'black'}, color=next(palette), xlim=(0,81), ylim=(0,1))#, ax=ax1)
	plt.xlabel('Average gaps (%)')
	plt.ylabel('Spearman correlation ({})'.format(r'$\rho$'))
	stats_patch = mpatches.Patch(color='white', label='$r^2$ = 0.001;  $pval$ = $1.1x10^-$$^1$$^8$')#, contains=False)
	plt.legend(handles=[stats_patch])
	plt.text(-3.4, 1.2, "d", fontsize=20, fontweight='bold', va='top', ha='right')
	plt.tight_layout()
	plt.show()

	next(palette)
	next(palette)
	df = df[df["theight_var"] <= 15]
	ax2 = sns.jointplot(x="theight_var", y=SCORES_LST[0], data=df, kind='reg', stat_func=pearsonr, line_kws={'color':'black'}, color=next(palette), xlim=(0,13), ylim=(0,1))#, ax=ax2)
	plt.xlabel('Deviation from ulrametricity')
	plt.ylabel("")
	stats_patch = mpatches.Patch(color='white', label='$r^2$ = 0.001;  $pval$ = $4.9x10^-$$^1$$^6$')  # , contains=False)
	plt.legend(handles=[stats_patch])
	plt.text(-0.65, 1.2, "e", fontsize=20, fontweight='bold', va='top', ha='right')
	plt.tight_layout()
	plt.show()






if __name__ == '__main__':
	dirpath = SUMMARY_FILES_DIR if platform.system() == 'Linux' else DATA_PATH
	#calc_empirical_features()


	df_val = pd.read_csv(dirpath + 'scores_per_ds_validation_with_more_atts.csv')
	df_train = pd.read_csv(dirpath + 'scores_per_ds_20_1_ytransformed_exp_with_more_atts.csv')
	#corr_plot(df_val)
	corr_plot_more_atts(df_train)
	#plot_distributions(df)
