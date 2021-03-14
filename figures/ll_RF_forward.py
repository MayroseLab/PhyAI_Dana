import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")

from defs import *

import itertools
sns.set_style("white")
fig, axarr = plt.subplots()
palette = itertools.cycle(sns.color_palette('colorblind'))




def plot_lines(lls_arr, RFs_arr):
	moveN = np.arange(0,16)
	'''
	fig, ax = plt.subplots()
	plt.plot(moveN, RFs_arr, '--r', label='Robinson-Foulds distance\n(from the maximum-likelihood tree)')
	plt.legend()
	#ax.tick_params('vals', colors='r')

	# Get second axis
	ax2 = ax.twinx()
	plt.plot(moveN, lls_arr, '--b', label='-Log-likelihood')
	plt.legend()
	#ax.tick_params('vals', colors='b')
	'''

	fig, ax1 = plt.subplots()

	color = next(palette)
	ax1.set_xlabel('Iteration number')
	ax1.set_ylabel('Robinson-Foulds distance\n(from the maximum-likelihood tree)', color=color)
	ax1.plot(moveN, RFs_arr, color=color)
	ax1.tick_params(axis='y', labelcolor=color)

	ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
	color = next(palette)
	ax2.set_ylabel('Log-likelihood', color=color)  # we already handled the x-label with ax1
	ax2.plot(moveN, lls_arr, color=color)
	ax2.tick_params(axis='y', labelcolor=color)

	plt.xticks(moveN)
	fig.tight_layout()
	plt.savefig("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/summary_files/" + "FigS5.tif", dpi=300)
	#plt.show()



if __name__ == '__main__':
	lls_arr = [-58580.26187, -58549.031018, -58510.397707, -58486.387171, -58465.672468, -58434.728162, -58422.113683, -58418.129932, -58414.237602,	-58410.68684, -58407.985743, -58407.342491, -58403.448113, -58403.447770000006, -58403.447758, -58403.447757]
	RFs_arr = [26, 22, 20, 18, 16, 14, 14, 10, 8, 8, 6, 4, 2, 2, 2, 0]
	plot_lines(lls_arr, RFs_arr)