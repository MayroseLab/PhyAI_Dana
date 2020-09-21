import numpy as np
from ete3 import *
import pandas as pd
#import regex as re
from collections import Counter
import numpy as np
from Bio import Phylo
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import *
from Bio.SeqRecord import SeqRecord
from io import StringIO
from collections import OrderedDict
#import matplotlib.pyplot as plt

PHYLIP_FORMAT = "phylip-relaxed"
SUBTREE1 = "subtree1"
SUBTREE2 = "subtree2"
TREES_PER_DS = "newicks_step1.csv"
LABEL = "d_ll_{}"
##################################################

def draw_tree(datapath):
	tree = "(Sp0003:1e-08,(Sp0020:0.0903209,((Sp0010:0.00975812,(Sp0009:0.0105802,Sp0016:0.00604999)N9:0.00301306)N6:0.00969131,(((Sp0024:0.00071452,(Sp0021:0.00139502,(Sp0007:0.0021314,Sp0027:0.00142406)N27:0.00071334)N19:0.00069954)N14:0.00717597,(Sp0000:0.0052456,(Sp0025:0.00191148,Sp0026:5e-08)N21:0.00761345)N15:0.00163044)N10:0.0167685,((Sp0011:0.458416,(Sp0005:0.0319714,(Sp0004:0.755483,Sp0022:0.0276043)N31:0.0249238)N23:0.0765557)N16:0.048788,((Sp0014:0.0111012,Sp0023:0.0153208)N24:0.00664935,(Sp0017:0.0128477,((Sp0006:0.00746376,(Sp0001:0.00180978,Sp0008:0.00540843)N43:0.00546766)N40:0.00786406,((Sp0002:0.0269194,(Sp0012:0.00356218,Sp0015:0.0093611)N49:0.00186054)N44:0.00414043,(Sp0018:0.00628908,Sp0013:0.0113868)N45:0.00226724)N41:0.00228271)N35:0.00467318)N25:0.00944122)N17:0.0206812)N11:0.0153946)N7:0.00336602)N5:0.0139865)N2:0.00239987,Sp0019:0.00122264);"
	tree_Sp2_N25 = "(Sp0003:1e-08,(Sp0020:0.0903209,((Sp0010:0.00975812,(Sp0009:0.0105802,Sp0016:0.00604999)N9:0.00301306)N6:0.00969131,(((Sp0024:0.00071452,(Sp0021:0.00139502,(Sp0007:0.0021314,Sp0027:0.00142406)N27:0.00071334)N19:0.00069954)N14:0.00717597,(Sp0000:0.0052456,(Sp0025:0.00191148,Sp0026:5e-08)N21:0.00761345)N15:0.00163044)N10:0.0167685,((Sp0011:0.458416,(Sp0005:0.0319714,(Sp0004:0.755483,Sp0022:0.0276043)N31:0.0249238)N23:0.0765557)N16:0.048788,((Sp0014:0.0111012,Sp0023:0.0153208)N24:0.00664935,(Sp0002:0.0269194,(Sp0017:0.0128477,((Sp0006:0.00746376,(Sp0001:0.00180978,Sp0008:0.00540843)N43:0.00546766)N40:0.00786406,((Sp0018:0.00628908,Sp0013:0.0113868)N45:0.00226724,(Sp0012:0.00356218,Sp0015:0.0093611)N49:0.00600097)N41:0.00228271)N35:0.00467318)N25:0.00472061)N44:0.00472061)N17:0.0206812)N11:0.0153946)N7:0.00336602)N5:0.0139865)N2:0.00239987,Sp0019:0.00122264);"
	tree_N44_N35 = "(Sp0003:1e-08,(Sp0020:0.0903209,((Sp0010:0.00975812,(Sp0009:0.0105802,Sp0016:0.00604999)N9:0.00301306)N6:0.00969131,(((Sp0024:0.00071452,(Sp0021:0.00139502,(Sp0007:0.0021314,Sp0027:0.00142406)N27:0.00071334)N19:0.00069954)N14:0.00717597,(Sp0000:0.0052456,(Sp0025:0.00191148,Sp0026:5e-08)N21:0.00761345)N15:0.00163044)N10:0.0167685,((Sp0011:0.458416,(Sp0005:0.0319714,(Sp0004:0.755483,Sp0022:0.0276043)N31:0.0249238)N23:0.0765557)N16:0.048788,((Sp0014:0.0111012,Sp0023:0.0153208)N24:0.00664935,(Sp0017:0.0128477,((Sp0002:0.0269194,(Sp0012:0.00356218,Sp0015:0.0093611)N49:0.00186054)N44:0.00414043,((Sp0006:0.00746376,(Sp0001:0.00180978,Sp0008:0.00540843)N43:0.00546766)N40:0.00786406,(Sp0018:0.00628908,Sp0013:0.0113868)N45:0.00454995)N35:0.00233659)N41:0.00233659)N25:0.00944122)N17:0.0206812)N11:0.0153946)N7:0.00336602)N5:0.0139865)N2:0.00239987,Sp0019:0.00122264);"
	
	d = pd.DataFrame.to_dict(pd.read_csv(datapath + "conversion_dict.csv"))['0']
	
	for k in d:
		kpad = "Sp" + str(k).zfill(4)
		orig_name = d[k].replace("_", " ")
		tree = tree.replace(kpad, orig_name)
		tree_Sp2_N25 = tree_Sp2_N25.replace(kpad, orig_name)
		tree_N44_N35 = tree_N44_N35.replace(kpad, orig_name)
	
	t, t2, t3 = Tree(tree, format=1), Tree(tree_Sp2_N25, format=1), Tree(tree_N44_N35, format=1)
	print(t.get_ascii(show_internal=True), "\n###########\n")
	print(t2.get_ascii(show_internal=True))
	
	ts = TreeStyle()
	ts.show_leaf_name = True
	# ts.rotation = 90
	ts.mode = "c"
	ts.arc_start = -180  # 0 degrees = 3 o'clock
	ts.arc_span = 180
	ts.scale = 200
	ts.branch_vertical_margin = 15
	ts.show_branch_length = True
	
	nst_general = NodeStyle()
	nst_general["fgcolor"] = "white"
	for tn in [t, t2, t3]:
		for l in tn.iter_leaves():
			l.img_style = nst_general
	
	nst1 = NodeStyle()
	# nst1["bgcolor"] = "LightSteelBlue"n #"Moccasin"
	nst1["fgcolor"] = "red"
	nst1["size"] = 5
	
	nst2 = NodeStyle()
	nst2["fgcolor"] = "blue"
	nst2["size"] = 5
	
	t_cp = t.copy()
	n2 = t & "Macaca nemestrina"
	n3 = t.get_common_ancestor("Chlorocebus aethiops", "Cercocebus atys", "Theropithecus gelada", "Papio anubis",
	                           "Macaca nemestrina", "Macaca mulatta", "Macaca fascicularis", "Mandrillus sphinx",
	                           "Mandrillus leucophaeus")
	n4 = t.get_common_ancestor("Macaca nemestrina", "Macaca mulatta", "Macaca fascicularis")
	n5 = t.get_common_ancestor("Cercocebus atys", "Theropithecus gelada", "Papio anubis", "Macaca nemestrina",
	                           "Macaca mulatta", "Macaca fascicularis", "Mandrillus sphinx", "Mandrillus leucophaeus")
	
	n2.set_style(nst1)
	n3.set_style(nst1)
	n4.set_style(nst2)
	n5.set_style(nst2)
	
	t.render(r"C:\Users\ItayMNB3\Dropbox\PhyloAI\\" + 'Fig3.pdf', tree_style=ts, dpi=3000)
	# t.render(r"C:\Users\ItayMNB3\Dropbox\PhyloAI\\" + 'fig_exp4_orig_Sp2_N25.PNG', tree_style=ts)
	# t2.render(r"C:\Users\ItayMNB3\Dropbox\PhyloAI\\" + 'fig_exp4_Sp2_N25.PNG', tree_style=ts)
	# t3.render(r"C:\Users\ItayMNB3\Dropbox\PhyloAI\\" + 'fig_exp4_N44_N35.PNG', tree_style=ts)
	# t_cp.render(r"C:\Users\ItayMNB3\Dropbox\PhyloAI\\" + 'fig_exp4_orig_N44_N35.PNG', tree_style=ts)





##########################

def prune_branch(t_orig, prune_name):
	'''
	get (a copy of) both subtrees after pruning
	'''
	t_cp_p = t_orig.copy()  				# the original tree is needed for each iteration
	prune_node_cp = t_cp_p & prune_name     # locate the node in the copied subtree
	assert prune_node_cp.up

	nname = prune_node_cp.up.name
	prune_loc = prune_node_cp
	prune_loc.detach()  # pruning: prune_node_cp is now the subtree we detached. t_cp_p is the one that was left behind
	t_cp_p.search_nodes(name=nname)[0].delete(preserve_branch_length=True)  # delete the specific node (without its childs) since after pruning this branch should not be divided
	if not nname:
		nname = "Nnew_p"

	return nname, prune_node_cp, t_cp_p


def regraft_branch(t_cp_p, rgft_node, prune_node_cp, rgft_name, nname):
	'''
	get a tree with the 2 concatenated subtrees
	'''
	# print("### rgft node: ", rgft_name)

	new_branch_length = rgft_node.dist /2
	t_temp = PhyloTree()  			   # for concatenation of both subtrees ahead, to avoid polytomy
	t_temp.add_child(prune_node_cp)
	t_curr = t_cp_p.copy()
	rgft_node_cp = t_curr & rgft_name  # locate the node in the copied subtree

	rgft_loc = rgft_node_cp.up
	rgft_node_cp.detach()
	t_temp.add_child(rgft_node_cp, dist=new_branch_length)
	t_temp.name = nname
	rgft_loc.add_child(t_temp, dist=new_branch_length)  # regrafting

	return t_curr

##################################################

if __name__ == '__main__':
	#datapath = r"D:\Users\Administrator\Dropbox\PhyloAI\data\training_datasets\\"
	datapath = r"C:\Users\ItayMNB3\Dropbox\PhyloAI\data\training_datasets\\"
	PHYML_TREE_FILENAME = "masked_species_real_msa.phy_phyml_tree_bionj.txt"
	MSA_PHYLIP_FILENAME = "masked_species_real_msa.phy"
	# felsenstein_log_likelihood(datapath + PHYML_TREE_FILENAME, datapath + MSA_PHYLIP_FILENAME, M, pai)
	# try_calc_leaves_features(tree, datapath + MSA_PHYLIP_FILENAME)
	#draw_tree(datapath)
	dirpath = os.getcwd()
	import raxml

	'''
	t_last_minus1_str = "(Sp0005:0.0299768,(Sp0011:0.60647,Sp0020:0.0851023)N2:0.0903137,((Sp0004:0.7,(((Sp0002:0.0219083,(Sp0017:0.0136272,(((Sp0018:0.00563617,Sp0013:0.0114737)N42:0.00224516,(Sp0006:0.00656193,(Sp0001:0.0013587,Sp0008:0.00528312)N51:0.00614451)N43:0.00574793)N34:0.00298905,(Sp0012:0.00523325,Sp0015:0.00616787)N35:0.00563218)N25:0.00346597)N17:0.00247621)N12:0.00719293,(Sp0014:0.0108875,Sp0023:0.0140792)N13:0.00740825)N10:0.0200846,(((Sp0010:0.00934953,(Sp0009:0.00929147,Sp0016:0.00579718)N27:0.00281884)N20:0.00963492,(Sp0003:1e-08,Sp0019:0.00115992)N21:0.0155618)N14:0.00272115,((Sp0024:0.00067441,(Sp0021:0.00132866,(Sp0007:0.00134969,Sp0027:0.0013507)N39:0.00067678)N31:0.00066422)N22:0.00683501,(Sp0000:0.00491566,(Sp0025:1e-08,Sp0026:5e-08)N33:0.00718346)N23:0.00147096)N15:0.0160162)N11:0.0151203)N9:0.702)N4:0.702,Sp0022:0.0468355)N3:6.099e-05);"
	t_last_minus1 = Tree(t_last_minus1_str, format=1)
	print(t_last_minus1.get_ascii(show_internal=True))
	
	t_str = "(Sp0005:0.0299768,(Sp0011:0.60647,(Sp0020:0.0851023,(((Sp0002:0.0219083,(Sp0017:0.0136272,(((Sp0018:0.00563617,Sp0013:0.0114737)N42:0.00224516,(Sp0006:0.00656193,(Sp0001:0.0013587,Sp0008:0.00528312)N51:0.00614451)N43:0.00574793)N34:0.00298905,(Sp0012:0.00523325,Sp0015:0.00616787)N35:0.00563218)N25:0.00346597)N17:0.00247621)N12:0.00719293,(Sp0014:0.0108875,Sp0023:0.0140792)N13:0.00740825)N10:0.0200846,(((Sp0010:0.00934953,(Sp0009:0.00929147,Sp0016:0.00579718)N27:0.00281884)N20:0.00963492,(Sp0003:1e-08,Sp0019:0.00115992)N21:0.0155618)N14:0.00272115,((Sp0024:0.00067441,(Sp0021:0.00132866,(Sp0007:0.00134969,Sp0027:0.0013507)N39:0.00067678)N31:0.00066422)N22:0.00683501,(Sp0000:0.00491566,(Sp0025:1e-08,Sp0026:5e-08)N33:0.00718346)N23:0.00147096)N15:0.0160162)N11:0.0151203)N9:1.15e-06)N5:0.0404572)N2:0.0903137,(Sp0004:1.4516,Sp0022:0.0468355)N3:6.099e-05);"
	t_last = PhyloTree(newick=t_str, alignment="C:\\Users\\ItayMNB3\\Dropbox\\PhyloAI\\" + MSA_PHYLIP_FILENAME,
	                   alg_format="iphylip", format=1)
	print(t_last.get_ascii(show_internal=True))
	'''
	'''
	# N10 , Sp0024
	t_str = "(Sp0005:0.0299768,(Sp0011:0.60647,(Sp0020:0.0851023,(((Sp0002:0.0219083,(Sp0017:0.0136272,(((Sp0018:0.00563617,Sp0013:0.0114737)N42:0.00224516,(Sp0006:0.00656193,(Sp0001:0.0013587,Sp0008:0.00528312)N51:0.00614451)N43:0.00574793)N34:0.00298905,(Sp0012:0.00523325,Sp0015:0.00616787)N35:0.00563218)N25:0.00346597)N17:0.00247621)N12:0.00719293,(Sp0014:0.0108875,Sp0023:0.0140792)N13:0.00740825)N10:0.0200846,(((Sp0010:0.00934953,(Sp0009:0.00929147,Sp0016:0.00579718)N27:0.00281884)N20:0.00963492,(Sp0003:1e-08,Sp0019:0.00115992)N21:0.0155618)N14:0.00272115,((Sp0024:0.00067441,(Sp0021:0.00132866,(Sp0007:0.00134969,Sp0027:0.0013507)N39:0.00067678)N31:0.00066422)N22:0.00683501,(Sp0000:0.00491566,(Sp0025:1e-08,Sp0026:5e-08)N33:0.00718346)N23:0.00147096)N15:0.0160162)N11:0.0151203)N9:1.15e-06)N5:0.0404572)N2:0.0903137,(Sp0004:1.4516,Sp0022:0.0468355)N3:6.099e-05);"
	#       "(Sp0005:0.0294646,(Sp0011:0.59499,(Sp0020:0.0934607,(((Sp0002:0.0217521,(Sp0017:0.0135506,(((Sp0018:0.00560708,Sp0013:0.0113141)N42:0.00226252,(Sp0006:0.00647705,(Sp0001:0.0013302,Sp0008:0.00525064)N51:0.00611572)N43:0.00581932)N34:0.00276965,(Sp0012:0.00516596,Sp0015:0.00614167)N35:0.00562778)N25:0.00343638)N17:0.00242913)N12:0.00704798,(Sp0014:0.0109123,Sp0023:0.0139947)N13:0.0074707)N10:0.0197984,(((Sp0010:0.00926167,(Sp0009:0.0092079,Sp0016:0.00573817)N27:0.00283838)N20:0.00945514,(Sp0003:1e-08,Sp0019:0.00115417)N21:0.0154632)N14:0.00277508,((Sp0024:0.00067077,(Sp0021:0.00132129,(Sp0007:0.00133827,Sp0027:0.00133801)N39:0.00067138)N31:0.00066003)N22:0.00681184,(Sp0000:0.00487686,(Sp0025:1e-08,Sp0026:5e-08)N33:0.00713307)N23:0.00145121)N15:0.0159483)N11:0.00740148)N9:0.00740148)N5:0.0353307)N2:0.0925657,(Sp0004:1.663,Sp0022:0.0201305)N3:0.0269788);"
	#t_last = Tree(t_str, format=1)
	t_last = PhyloTree(newick=t_str, alignment="C:\\Users\\ItayMNB3\\Dropbox\\PhyloAI\\" + MSA_PHYLIP_FILENAME, alg_format="iphylip", format=1)
	print(t_last.get_ascii(show_internal=True))
	#t_last_minus1_str = "(Sp0005:0.0299768,(Sp0011:0.60647,Sp0020:0.0851023)N2:0.0903137,((Sp0004:0.7,(((Sp0002:0.0219083,(Sp0017:0.0136272,(((Sp0018:0.00563617,Sp0013:0.0114737)N42:0.00224516,(Sp0006:0.00656193,(Sp0001:0.0013587,Sp0008:0.00528312)N51:0.00614451)N43:0.00574793)N34:0.00298905,(Sp0012:0.00523325,Sp0015:0.00616787)N35:0.00563218)N25:0.00346597)N17:0.00247621)N12:0.00719293,(Sp0014:0.0108875,Sp0023:0.0140792)N13:0.00740825)N10:0.0200846,(((Sp0010:0.00934953,(Sp0009:0.00929147,Sp0016:0.00579718)N27:0.00281884)N20:0.00963492,(Sp0003:1e-08,Sp0019:0.00115992)N21:0.0155618)N14:0.00272115,((Sp0024:0.00067441,(Sp0021:0.00132866,(Sp0007:0.00134969,Sp0027:0.0013507)N39:0.00067678)N31:0.00066422)N22:0.00683501,(Sp0000:0.00491566,(Sp0025:1e-08,Sp0026:5e-08)N33:0.00718346)N23:0.00147096)N15:0.0160162)N11:0.0151203)N9:0.702)N4:0.702,Sp0022:0.0468355)N3:6.099e-05);"
	#t_last_minus1 = Tree(t_last_minus1_str, format=1)
	#print(t_last_minus1.get_ascii(show_internal=True))
	nname, subtree1, subtree2 = prune_branch(t_last, "N10")
	rgft_node = subtree2 & "Sp0024"
	full_tree = regraft_branch(subtree2, rgft_node, subtree1, "Sp0024", nname)
	print(full_tree.get_ascii(show_internal=True))
	full_tree.write(format=1, outfile="C:\\Users\\ItayMNB3\\Dropbox\\PhyloAI\\" + "minus1_tree.txt")
	'''
	
	'''
	from sklearn import preprocessing
	import seaborn as sns
	palette = itertools.cycle(sns.color_palette('colorblind'))
	
	#df = pd.read_excel(datapath + "example4_preds.xlsx")
	df = pd.read_csv(r"C:\\Users\ItayMNB3\Dropbox\PhyloAI\\" + "with_preds_merged_28_1_example4b_ytransformed_exp10.csv")
	#df = pd.read_excel(r"C:\\Users\ItayMNB3\Dropbox\PhyloAI\example4_preds")
	true_dll = df["d_ll_prune"]
	pred_dll = df["pred"]  #df["pre_inverse_trans"]
	
	
	#scaler = preprocessing.StandardScaler()
	#scaled_x1 = scaler.fit_transform(true_dll.values.reshape(-1, 1))
	#scaled_x2 = scaler.fit_transform(pred_dll.values.reshape(-1, 1))
	#df["d_ll_prune"] = scaled_x1
	#df["pre_inverse_trans"] = scaled_x2
	
	ax = sns.distplot(true_dll.values, kde=False, label="big",color=next(palette), hist_kws=dict(alpha=0.5))#, bins=10)
	sns.distplot(pred_dll.values, kde=False, label="big",color=next(palette), hist_kws=dict(alpha=0.5))# , bins=10)
	
	plt.show()
	'''