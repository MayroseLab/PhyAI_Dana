import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")
import warnings
warnings.filterwarnings("ignore")			# TEMP

from defs import *

import seaborn as sns
sns.set_style("white")
horlines = [0, 1.16, 2.3, 4.6]


def plot_violin(d1, d2, horizontal=False):
	fig, axarr = plt.subplots(1, 2, sharey=True)
	if not horizontal:
		#fig, axarr = plt.subplots(1, 1, sharey=True)
		#ax = sns.violinplot(list(d1.keys())*len(list(d1.values())[0]), list(d1.values())[0], inner="point", scale="width", ax=axarr[0]) #, color="blue")
		sns.violinplot(list(d1.keys()) * len(list(d1.values())[0]), list(d1.values())[0], inner="quartiles", scale="width",ax=axarr[0]) #, orient='h'
		sns.violinplot(list(d2.keys())*len(list(d2.values())[0]), list(d2.values())[0], inner="quartiles", scale="width", ax=axarr[1], color="darkred")

		#fig.suptitle("Random-Forest performance evaluation", y=1, size=16)
		axarr[0].set_ylabel("rank", size=14)
		#ax.set_ylabel("rank", size=14)
		plt.ylim([-5,117])
		plt.yticks(np.array([1,10,20,30,40,50,60,70,80,90,100,110,117]))

		fig.tight_layout()
		plt.savefig(DATA_PATH + 'violinplot.png')
		plt.show()

	else:
		sns.violinplot(list(d1.values())[0], list(d1.keys()) * len(list(d1.values())[0]), inner="quartiles",
					   scale="width", ax=axarr[0])  # , orient='h'
		sns.violinplot(list(d2.values())[0],list(d2.keys()) * len(list(d2.values())[0]), inner="quartiles",
					   scale="width", ax=axarr[1], color="darkred")

		axarr[0].set_xlabel("rank", size=14)
		plt.xlim([-5, 117])
		plt.xticks(np.array([1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 117]))

		fig.tight_layout()
		plt.savefig(DATA_PATH + 'violinplot.png')
		plt.show()


def plot_hist(d1, d2):
	fig, axarr = plt.subplots(1, 1)
	sns.distplot(list(d1.values())[0], kde=False, color='darkred')
	#sns.kdeplot(list(d1.values())[0], shade=True, color='darkred')
	plt.xticks(np.array([1,2,3,4,5, 10, 20, 30]))
	plt.yticks(np.arange(0,1001,200))
	plt.xlim([0.7, 30])


	axarr.set_xlabel("Rank", size=18)
	axarr.set_ylabel("% datasets", size=18)

	fig.tight_layout()
	#plt.savefig(DATA_PATH + 'hist.png')
	plt.show()



if __name__ == '__main__':
	# test
	pred_in_true = [1, 6, 4, 4, 6, 1, 4, 5, 3, 4, 2, 6, 48, 7, 1, 12, 2, 4, 4, 7, 2, 6, 1, 1, 2, 1, 1, 1, 2, 6, 1, 3, 2, 3, 2, 4, 3, 6, 2, 3, 27, 1, 1, 1, 4, 6, 27,
	                1, 10, 4, 6, 1, 3, 23, 13, 5, 2, 1, 3, 3, 3, 12, 2, 5, 2, 2, 3, 2, 5, 3, 15, 2, 2, 4, 5, 5, 2, 7, 2, 4, 4, 2, 3, 6, 1, 1, 3, 2, 7, 5]
	true_in_pred = [1, 8, 16, 31, 17, 1, 5, 23, 8, 7, 3, 5, 4, 7, 1, 3, 4, 4, 14, 19, 6, 19, 1, 1, 3, 1, 1, 1, 3, 49, 1, 5, 3, 3, 7, 3, 12, 18, 3, 2, 5, 1, 1, 1,
	                3, 21, 23, 1, 7, 28, 21, 1,10, 28, 5, 3, 17, 1, 31, 3, 7, 2, 7, 2, 15,3, 19, 7, 19, 19, 4, 9, 2, 6, 3, 3, 2, 3, 4, 3,10, 2, 2, 7, 1, 1, 2, 12, 19, 7]

	#plot_violin({"rank of best predicted\nSPR move": pred_in_true}, {"predicted rank of best\nSPR move":true_in_pred})
	plot_hist(pred_in_true, true_in_pred)