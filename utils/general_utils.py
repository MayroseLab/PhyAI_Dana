import csv
import sys, logging
from Bio import AlignIO, Alphabet
from defs import *
from utils.msa_functions import get_msa_properties


def change_path_permissions_to_777(path):
	os.chmod(path, 0o777)
	for root, dirs, files in os.walk(path):
		for dir in dirs:
			try:
				os.chmod(os.path.join(root, dir), 0o777)
			except:
				pass
		for file in files:
			try:
				os.chmod(os.path.join(root, file), 0o777)
			except:
				pass



def convert_fasta_to_phylip(input_file, output_file):
	with open(input_file, "rU") as input_handle:
		alignments = AlignIO.parse(input_handle, FASTA_FORMAT)
		with open(output_file, "w") as output_handle:
			AlignIO.write(alignments, output_handle, PHYLIP_FORMAT)

def convert_phylip_to_nexus(input_file, output_file):
	AlignIO.convert(input_file, PHYLIP_FORMAT, output_file, NEXUS_FORMAT, alphabet=Alphabet.generic_dna)


def convert_phylipInterleaved_to_sequential_relaxed(msa_file, output_file):
	with open(msa_file, "rU") as input_handle:
		alignments = AlignIO.read(input_handle, PHYLIP_FORMAT)

	with open(output_file, "w") as output_handle:
		out_handle = AlignIO.PhylipIO.SequentialPhylipWriter(output_handle)
		out_handle.write_alignment(alignments, id_width=30)


def init_commandline_logger(logger):
	logger.setLevel(logging.DEBUG)
	# create console handler and set level to debug
	ch = logging.StreamHandler(sys.stdout)
	ch.setLevel(logging.DEBUG)
	# create formatter
	formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	# add formatter to ch
	ch.setFormatter(formatter)
	# add ch to logger
	logger.addHandler(ch)


def num_of_params_in_model(model):
	cnt = 0
	if model.endswith("+I+G"):
		cnt += 2
	elif model.endswith("G") or model.endswith("I"):
		cnt += 1

	if model.startswith("JC") or model.startswith("F81"):
		cnt += 1
	elif model.startswith("K2P") or model.startswith("HKY"):
		cnt += 2
	elif model.startswith("GTR") or model.startswith("SYM"):
		cnt += 6

	if model.startswith("JC") or model.startswith("K2P") or model.startswith("SYM"):
		cnt += 1
	elif model.startswith("F81") or model.startswith("HKY") or model.startswith("GTR"):
		cnt += 3
	return cnt


def is_aicc(ntaxa, nchars, model):
	#print(ntaxa, nchars, model)
	nparams = 2*ntaxa-3 + num_of_params_in_model(model)
	return nchars/nparams < 40


def split_model_components(model):
	pinvar = gamma = False
	if model.endswith("+I+G"):
		pinvar = True
		gamma = True
	elif model.endswith("+G"):
		gamma = True
	elif model.endswith("+I"):
		pinvar = True

	base_model = re.sub("\+.*$", "", model)
	return base_model, pinvar, gamma





def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		return False


def reorderLegend(ax=None,order=None,unique=False):
	"""
	Returns tuple of handles, labels for axis ax, after reordering them to conform to the label order `order`, and if unique is True, after removing entries with duplicate labels.
	"""
	if ax is None: ax=plt.gca()

	handles, labels = ax.get_legend_handles_labels()
	# sort both labels and handles by labels
	labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
	if order is not None: # Sort according to a given list (not necessarily complete)
		keys=dict(zip(order,range(len(order))))
		labels, handles = zip(*sorted(zip(labels, handles), key=lambda t,keys=keys: keys.get(t[0],np.inf)))
	if unique: # Keep only the first of each handle
		labels, handles= zip(*unique_everseen(zip(labels,handles), key = labels))
	ax.legend(handles, labels)

	return (handles, labels)


def unique_everseen(seq, key=None):
	seen = set()
	seen_add = seen.add
	return [x for x,k in zip(seq,key) if not (k in seen or seen_add(k))]


def is_file_empty(filepath):
	"""
	:param filepath:
	:return: True if filepath doesn't exist or is empty
	"""
	if os.path.exists(filepath):
		with open(filepath) as fpr:
			if re.search("\S", fpr.read()):
				return False
	return True