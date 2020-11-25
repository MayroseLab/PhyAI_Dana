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


def calc_empirical_features():
	#df = pd.read_csv(dirpath + CHOSEN_DATASETS_FILENAME).reset_index(drop=True)
	#df2 = pd.read_csv(dirpath + SCORES_PER_DS.format("20_1_4200_ytransformed_exp")).reset_index(drop=True)
	df = pd.read_csv('/groups/itay_mayrose/danaazouri/PhyAI/validation_set2/summary_files/sampled_datasets.csv').reset_index(drop=True)
	df2 = pd.read_csv("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/summary_files/results_saturation/scores_per_ds_20_1_validation_set_4200_ytransformed_exp_orig.csv").reset_index(drop=True)
	for index, row in df2.iterrows():
		path = row["path"]
		sum_ds_df = SUMMARY_PER_DS.format(path, "prune", "br", "1")
		if os.path.exists(sum_ds_df):
			df2.loc[index, "ntaxa"] = df.loc[df["path"] == path, "ntaxa"].values[0]
			df2.loc[index, "nchars"] = df.loc[df["path"] == path, "nchars"].values[0]
			df2.loc[index, "tbl"] = get_total_branch_lengths(path + PHYML_TREE_FILENAME.format('bionj'))
			df2.loc[index, "ll"] = pd.read_csv(sum_ds_df).loc[0, "orig_ds_ll"]
		
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
	#fig, (ax1, ax2, ax3) = plt.subplots(1, 3)

	ax1 = sns.jointplot(x="ntaxa", y=SCORES_LST[0], data=df, kind='reg', stat_func=pearsonr, line_kws={'color':'black'}, color=next(palette),xlim=(5, 70.3), ylim=(0,1))#, ax=ax1)
	plt.xlabel('Number of taxa')
	plt.ylabel('Spearman correlation ({})'.format(r'$\rho$'))
	#stats_patch = mpatches.Patch(color='white', label='$r^2$ = 0.0002;  $pval$ = 0.34')#, contains=False)
	stats_patch = mpatches.Patch(color='white', label='$r^2$ = 0.009;  $pval$ = $7.8x10^-$$^1$$^0$')#, contains=False)
	plt.legend(handles=[stats_patch])
	plt.text(0.11, 1.35, "a", fontsize=20, fontweight='bold', va='top', ha='right')
	plt.tight_layout()
	plt.show()

	df = df[df["tbl"] < 30]
	ax2 = sns.jointplot(x="tbl", y=SCORES_LST[0], data=df, kind='reg', stat_func=pearsonr, line_kws={'color':'black'}, color=next(palette),xlim=(-0.3, 31), ylim=(0,1))#, ax=ax2)

	plt.xlabel('Total branch lengths')
	plt.ylabel("")
	# plt.xscale("log")
	#stats_patch = mpatches.Patch(color='white', label='$r^2$ = 0.007;  $pval$ = $1.1x10^-$$^8$')#, contains=False)
	stats_patch = mpatches.Patch(color='white', label='$r^2$ = 0.025;  $pval$ = $3.1x10^-$$^2$$^4$')  # , contains=False)
	plt.legend(handles=[stats_patch])
	plt.text(-0.9, 1.35, "b", fontsize=20, fontweight='bold', va='top', ha='right')
	plt.tight_layout()
	plt.show()
	
	df_val = pd.read_csv(dirpath + SCORES_PER_DS.format("20_1_validation_set_4200_ytransformed_exp"))
	#df_val.replace([df_val["Database"] == "orthoMam"]["Database"], "OrthoMaM")
	df_val["Database"] = df_val["Database"].replace(["orthoMam"], "OrthoMaM")
	df_val["Set"] = "validation"
	df["Set"] = "training"
	df_inc_val = pd.concat([df, df_val])
	ax3 = sns.boxplot(x="Database", y=SCORES_LST[0], data=df_inc_val, showfliers=False,order=["ploiDB","protDBs", "selectome", "TreeBASE", "OrthoMaM", "PANDIT"], color="lightgrey")#, ax=ax3)
	plt.ylim(0,1)
	plt.ylabel("")
	plt.setp(ax3.get_xticklabels(), rotation=20)
	plt.text(-0.8, 1.09, "c", fontsize=20, fontweight='bold', va='top', ha='right')
	plt.tight_layout()
	plt.show()

	gpd_df = df_inc_val.groupby('Database')
	for i, df_db in gpd_df:
		print(df_db['Database'].values[0])
		print(df_db[SCORES_LST[0]].values.mean())
	
	
	
	#plt.tight_layout()
	#plt.show()
	
	
	
	'''
	gs0 = gridspec.GridSpec(2, 2)
	gs00 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[0, :])
	gs01 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[1, :])
	
	ax1, ax2, ax3 = f.add_subplot(gs00[:]), f.add_subplot(gs01[:, :1]), f.add_subplot(gs01[:, 1:2])
	
	df_inc_val = pd.concat([df, pd.read_csv(dirpath + SCORES_PER_DS.format("26_1_validation_set_ytransformed_exp"))])

	#sns.catplot(x="Database", y=SCORES_LST[0], data=df, ax=ax1, order=[DBS[0],DBS[2],DBS[3],DBS[1]], kind='swarm)
	sns.boxplot(x="Database", y=SCORES_LST[0], data=df_inc_val, ax=ax1, showfliers=False, order=["protDBs", "PANDIT", "ploiDB", "orthoMam", "TreeBASE", "selectome"])
	ax1.legend(loc='best', frameon=True)
	
	sns.jointplot(x='ntaxa', y =SCORES_LST[0], data = df, kind ='reg', stat_func=pearsonr, color=next(palette), ax=ax2)
	#sns.regplot(x="ntaxa", y=SCORES_LST[0], data=df, ax=ax2, color=next(palette))
	sns.regplot(x="tbl", y=SCORES_LST[0], data=df, ax=ax3, color=next(palette))
	ax3.set_xlim(0,10)
	
	#ax1.set(xlabel='Spearman correlation ({})'.format(r'$\rho$'), ylabel="% empirical datasets", title='score 1',
	#       yticklabels=yperc1, xlim=(0, 1))
	
	f.tight_layout()
	plt.show()
	'''
	



if __name__ == '__main__':
	dirpath = SUMMARY_FILES_DIR if platform.system() == 'Linux' else DATA_PATH
	#calc_empirical_features()
	
	
	df_val = pd.read_csv(dirpath + 'scores_per_ds_validation_with_atts.csv')
	df_train = pd.read_csv(dirpath + 'scores_per_ds_20_1_ytransformed_exp_with_atts.csv')
	#plot_distributions(df)
	corr_plot(df_train)
	#binned_att_to_scores(df)