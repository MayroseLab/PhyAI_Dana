import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")

from defs import *



def traverse_data_dirs(func, csv_source=SUMMARY_FILES_DIR+CHOSEN_DATASETS_FILENAME, slice=(False,False), *args):
	"""
	:param func: the first argument of func should be the cluster path
	:param args:
	:return:
	"""

	func_outputs = []

	df = pd.read_csv(csv_source)
	if slice[0]:
		start = slice[0]
		df = df.loc[int(start):(int(start)+int(slice[1]) -1)]
	for index, row in df.iterrows():
		func_outputs.append(func(row["path"], *args))


	return func_outputs