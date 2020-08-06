import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")
import warnings
warnings.filterwarnings("ignore")			# TEMP

from defs import *
from ete3 import Tree, TreeStyle
import csv
from Bio.Blast import NCBIWWW



def layout(node):
	sp_dict_path = DATA_PATH + "ids_conversion_map.csv"
	with open(sp_dict_path) as csv_file:
		reader = csv.reader(csv_file)
		mydict = dict(reader)
		code2desc = {v: k for k, v in mydict.items()}

	# If node is a leaf, add the nodes name and a its scientific name
	if node.is_leaf():
		descFace = faces.TextFace(code2desc[node.name], fsize=10)
		faces.add_face_to_node(descFace, node, column=0) #, aligned=True)


def get_species_name(node_name_string):
	sp_dict_path = DATA_PATH + "ids_conversion_map.csv"
	with open(sp_dict_path) as csv_file:
		reader = csv.reader(csv_file)
		mydict = dict(reader)
		code2name = {v: k for k, v in mydict.items()}
	return code2name[node_name_string]



def draw_ranked_edges(df, t):
	ts = TreeStyle()
	ts.show_leaf_name = False
	ts.rotation = 90
	ts.mode = "c"
	ts.arc_start = -180  # 0 degrees = 3 o'clock
	ts.arc_span = 180
	ts.scale = 500
	ts.branch_vertical_margin = 7

	nstyle = NodeStyle()
	nstyle["shape"] = "sphere"
	nstyle["size"] = 1
	nstyle["fgcolor"] = "black"

	for i,n in enumerate(t.traverse()):
		n.set_style(nstyle)
		if i == 0:
			continue
		true = int(df.loc[df["prune_name"] == n.name, "ranked_label"].values[0])
		pred = int(df.loc[df["prune_name"] == n.name, "ranked_pred"].values[0])
		T = TextFace(str(true))
		P = TextFace(str(pred))
		T.background.color = "LightGreen"
		P.background.color = "orange"
		T.margin_left, T.margin_right  = 10, 10
		P.margin_left, P.margin_right = 10, 10
		n.add_face(T, column=0, position="branch-top")
		n.add_face(P, column=1, position="branch-bottom")

	#ts.layout_fn = layout
	t.show(tree_style=ts)
	t.render(DATA_PATH + 'fig2.PNG', tree_style=ts)





if __name__ == '__main__':
	df = pd.read_csv(DATA_PATH + "1.csv")
	#tree_path = DATA_PATH + PHYML_TREE_FILENAME.format("bionj")
	tree_path = DATA_PATH + "tree.txt"
	msa_path = DATA_PATH + MSA_PHYLIP_FILENAME

	t = Tree(newick=tree_path, format=1)
	#t = PhyloTree(newick=tree_path, alignment=msa_path, alg_format="iphylip", format=1)
	#t.set_species_naming_function(get_species_name)
	draw_ranked_edges(df, t)


