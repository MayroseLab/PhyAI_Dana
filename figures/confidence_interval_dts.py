import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")
import warnings
warnings.filterwarnings("ignore")			# TEMP

from defs import *




def plot_pred_true(X, all_preds, all_true, errs_down, errs_up):
	fig, ax = plt.subplots(1)

	#print(len(all_preds), len(all_true), len(errs_down), len(errs_up))
	ax.plot(X, all_true, 'bo', label='true values')
	#ax.plot(X, all_preds, 'ro', label='predicted values')
	#plt.fill_between(X, errs_down, errs_up, color='gray', alpha=0.2)
	ax.errorbar(x=X, y=all_preds, yerr=[errs_down,errs_up], label='predicted values',  fmt='ro')
	plt.xlabel('dataset no.')
	plt.legend(loc='upper right')
	plt.show()