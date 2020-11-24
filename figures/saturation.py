import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")
import warnings
warnings.filterwarnings("ignore")			# TEMP
import itertools
import seaborn as sns

from defs import *
import scipy.stats as stats
N_DATASETS_COL = "n_datasets"

sns.set_style("white")
fig, axarr = plt.subplots()
#palette = itertools.cycle(sns.color_palette('pastel'))




def plot_scores(df):
	#df = df.set_index('ndata')
	plt.figure()
	ax = plt.gca()
	#palette = itertools.cycle(sns.color_palette('colorblind'))
	leg = ax.legend()

	df_grouped = df.groupby('move_type', sort=True)
	for j, (name, group) in enumerate(df_grouped):
		for rank_mode in [False,True]:
			rank_label = "" if not rank_mode else "_ranked"
			subdf = group[group['if_ranked'] == rank_mode]
			color = next(palette)
			subdf.plot(x='ndata', y='2scores_mean', color=color, label=name+rank_label, ax=ax)
			subdf.plot.scatter(x='ndata', y='2scores_mean', color=color, ax=ax)
			subdf.plot(x='ndata', y='percentages', color=color, label=name + rank_label + "_precentages", linestyle='--', ax=ax)
			subdf.plot.scatter(x='ndata', y='percentages', color=color, ax=ax)

	ax.legend(loc='right', frameon=True)
	plt.xticks([500, 1000, 1500, 2000])
	plt.ylabel('score\n(mean of both scores)')
	#plt.title('Here data were sampled randomly, ntaxa 15 and 30')
	#plt.savefig(DATA_PATH + 'saturation_no_60_per.png')
	
	plt.show()



def plot_scores2(df):
	#palette = itertools.cycle(sns.color_palette('colorblind'))
	axarrs_locs = [(0, 0), (0, 1), (1, 0), (1, 1)]
	ylims = [(0.5, 1), (1, 40), (1, 40), (0, 30)]

	for i, loc in enumerate(axarrs_locs):
		color = next(palette)
		ax = axarr[axarrs_locs[i]]
		ax.set_ylim(ylims[i])
		sns.boxplot(x=N_DATASETS_COL, y=SCORES_LST[i], data=df, ax=ax, color=color, showfliers=False)
	
	fig.tight_layout()
	plt.show()
	#plt.savefig(dirpath + "saturation_plot.png")


def plot_corr(df):
	sns.set_context("paper", font_scale=1.3)
	corr_col = SCORES_LST[0]
	#df[corr_col] = df[corr_col].apply(np.sqrt)
	ax = sns.boxplot(x=N_DATASETS_COL, y=corr_col, data=df, showfliers=False, palette="husl", saturation=0.6, linewidth=0.6)
	ax.set_ylim(0,1)

	plt.ylabel('Spearman correlation ({})'.format(r'$\rho$'), size = 12)
	plt.xlabel('Number of empirical datasets in training', size = 12)
	fig.tight_layout()
	fig.set_size_inches(7, 4, forward=True)
	#plt.savefig("C:\\Users\\ItayMNB3\\Dropbox\\PhyloAI\\PhyAI_writing\\figures\\" + "FigS3.tif", dpi=300)
	plt.savefig(dirpath + "saturation.tif", dpi=300)
	plt.show()


	#plt.savefig(dirpath + "saturation_plot.png")


def concat_n_features(dirpath, ndots, xticks):
	'''
	pre-merge iterations results to one csv with "n_datasets" col
	'''
	df = pd.read_csv(dirpath + SCORES_PER_DS.format("20_1_{}_ytransformed_exp".format(ndots[0])))
	df[N_DATASETS_COL] = xticks[0]
	for i,n in enumerate(ndots[1:]):
		csv = dirpath + SCORES_PER_DS.format("20_1_{}_ytransformed_exp".format(n))
		df2 = pd.read_csv(csv)
		df2[N_DATASETS_COL] = xticks[i+1]
		df = pd.concat([df, df2])
	
	return df





if __name__ == '__main__':
	dirpath = SUMMARY_FILES_DIR if platform.system() == 'Linux' else DATA_PATH
	#df = pd.read_csv(dirpath + "results_summary.csv")
	#df = pd.read_csv(dirpath + "results_summary.csv")   # _by_ntaxa   , _sr
	#plot_scores
	
	ndots_older = [600, 840, 1310, 1790, 2260, 2670, 3210, 3870, 3860, 4670, 6060]
	xticks_older = [600, 900, 1300, 1800, 2300,2700, 3200, 3700, 4000, 4700, 6000]
	
	ndots = ['1500_k2', 2200, 3000, 4000, 5000, 5850]
	xticks = [1500, 2100, 3000, 4200, 5000, 6000]
	dirpath += 'results_saturation/'
	df = concat_n_features(dirpath, ndots, xticks)
	#plot_scores2(df)
	plot_corr(df)


	a1500 = df.loc[df[N_DATASETS_COL] == xticks[0], SCORES_LST[0]].dropna().values
	a2200 = df.loc[df[N_DATASETS_COL] == xticks[1], SCORES_LST[0]].dropna().values
	a3000 = df.loc[df[N_DATASETS_COL] == xticks[2], SCORES_LST[0]].dropna().values
	a4000 = df.loc[df[N_DATASETS_COL] == xticks[3], SCORES_LST[0]].dropna().values
	a5000 = df.loc[df[N_DATASETS_COL] == xticks[4], SCORES_LST[0]].dropna().values
	a6000 = df.loc[df[N_DATASETS_COL] == xticks[5], SCORES_LST[0]].dropna().values

	F, p = stats.f_oneway(a3000, a6000)
	print(p)