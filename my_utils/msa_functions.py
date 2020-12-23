import os
from Bio import AlignIO
from defs import *
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment



def get_msa_from_file(msa_file_path):
	#open file if exists
	if not os.path.exists(msa_file_path):
		return None
	try:
		msa = AlignIO.read(msa_file_path, PHYLIP_FORMAT)
	except:
		return None
	return msa


def get_msa_properties(msa):
	"""
	:param msa: bio.AlignIO format or path to msa file
	:return:
	"""
	if isinstance(msa, str):
		msa = get_msa_from_file(msa)
	ntaxa = len(msa)
	nchars = msa.get_alignment_length()

	return ntaxa, nchars


def remove_masked_sites(msa, dest_filename=None):
	"""
	:param msa: bio.AlignIO format or path to msa file
	:return:
	"""
	if isinstance(msa, str):
		msa = get_msa_from_file(msa)
	ntaxa = len(msa)
	nchars = msa.get_alignment_length()
	new_msa = copy.deepcopy(msa)
	for col_i in range(nchars - 1, -1, -1):
		col = msa[:, col_i]
		if not re.search("[ACGTacgt]", col, re.IGNORECASE):
			new_msa = new_msa[:, :col_i] + new_msa[:, col_i + 1:]
	if dest_filename:
		with open(dest_filename, "w") as fpw:
			AlignIO.write(new_msa, fpw, format=PHYLIP_FORMAT)
	return new_msa


def remove_gapped_sites(msa, thres):
	"""
	removes sites that contain more than thres gaps in a column
	:param msa: a Bio.Align.MultipleSeqAlignment object
	:param thres: a floating number in [0,1], the maximal fraction of
	sites that may hold gaps within every alignment column
	:return: a new msa without sites that have more than thres gaps
	"""
	ntaxa = len(msa)
	nchars = msa.get_alignment_length()
	new_msa = copy.deepcopy(msa)
	for col_i in range(nchars - 1, -1, -1):
		if msa[:, col_i].count("-") / ntaxa > thres:
			new_msa = new_msa[:, :col_i] + new_msa[:, col_i + 1:]
	return new_msa

def rewrite_in_phylip(msa_file):
	data = []
	with open(msa_file, 'r') as fp:
		for line in fp.readlines():
			nline = line.rstrip("\r\n")
			re_name_only = re.search("^(\S+)\s+\S+",nline)
			if re_name_only:
				name_only = re_name_only.group(1)
				end_name_ix = len(name_only) +1
				with_spaces = nline[:end_name_ix] + "      " + nline[end_name_ix:]
				nline = with_spaces
			data.append(nline)

	with open(msa_file, 'w') as nf:
		nf.write("\n".join(data))


def get_seqs_dict(msa_file):
	alignment = get_msa_from_file(msa_file)
	seqs_dict = {seq.id: seq.seq for seq in alignment}
	return seqs_dict


def trunc_msa2(subt, alignment, trunc_msa_path=None):

	if type(subt) == str:
		subt = Tree(newick=subt, format=1)
	names = subt.get_leaf_names()
	new_msa = []
	all_names = [rec.name for rec in list(alignment)]
	for name in names:
		new_msa.append(alignment[all_names.index(name), :])
	
	trunc_msa = MultipleSeqAlignment(new_msa)
	
	if trunc_msa_path:
		AlignIO.write(trunc_msa, trunc_msa_path, PHYLIP_FORMAT)
		rewrite_in_phylip(trunc_msa_path)

	return trunc_msa


def trunc_msa(subt, seqs_dict, trunc_msa_path=None):
	if type(subt) == str:
		subt = Tree(newick=subt, format=1)
	records = []
	for leaf in subt.iter_leaves():  # go over all leaves (all species) in this subtree
		leaf_name = leaf.name
		records.append(SeqRecord(seqs_dict[leaf_name], id=leaf_name))
	trunc_msa = MultipleSeqAlignment(records)
	if trunc_msa_path:
		AlignIO.write(trunc_msa, trunc_msa_path, PHYLIP_FORMAT)
		rewrite_in_phylip(trunc_msa_path)

	return trunc_msa


def hamming_distance(seq1, seq2):
	cnt = 0
	for i in range(len(seq1)):
		cnt += int(seq1[i]!=seq2[i])
	return cnt



def mask_species_names_in_msa(msa_file, dest_msa_file, dest_dir):
	msa = get_msa_from_file(msa_file)

	conversion_dict = {}
	for i, seq in enumerate(msa):
		orig_id = seq.id
		new_id = "Sp" + str(i).zfill(4)
		conversion_dict[new_id] = orig_id
		seq.id = new_id

	with open(dest_msa_file, "w") as output_handle:
		AlignIO.write(msa, output_handle, PHYLIP_FORMAT)

	df = pd.DataFrame.from_dict(conversion_dict, orient='index')
	df.to_csv(dest_dir + "conversion_dict.csv")
	return





if __name__ == '__main__':
	#print(remove_masked_sites("real_msa.phy"))
	#dirpath = r"D:\Users\Administrator\Dropbox\PhyloAI\data\training_datasets\\"
	#mask_species_names_in_msa(dirpath + "1.Rhopalaea.phy", dirpath + "1.xx.phy", dirpath)
	
	alignment = get_msa_from_file(DATA_PATH+MSA_PHYLIP_FILENAME)
	trunc_msa2(Tree(DATA_PATH + PHYML_TREE_FILENAME.format('bionj'), format=1), alignment, DATA_PATH + "test.phy")