import re, os, shutil
import pandas as pd
from ete3 import Tree
from subprocess import Popen, PIPE, STDOUT

RAXML_NG_SCRIPT = "raxml-ng"  # after you install raxml-ng on your machine
MSA_PHYLIP_FILENAME = "masked_species_real_msa.phy"


def prune_branch(t_orig, prune_name):
    '''
    returns (a copy of) both ETE subtrees after pruning
    '''
    t_cp_p = t_orig.copy()  # the original tree is needed for each iteration
    assert t_cp_p & prune_name  # todo Oz: add indicative error
    prune_node_cp = t_cp_p & prune_name  # locate the node in the copied subtree
    assert prune_node_cp.up
    
    nname = prune_node_cp.up.name
    prune_loc = prune_node_cp
    prune_loc.detach()  # pruning: prune_node_cp is now the subtree we detached. t_cp_p is the one that was left behind
    t_cp_p.search_nodes(name=nname)[0].delete(
        preserve_branch_length=True)  # delete the specific node (without its childs) since after pruning this branch should not be divided
    return nname, prune_node_cp, t_cp_p


def regraft_branch(t_cp_p, prune_node_cp, rgft_name, nname, preserve=False):
    '''
    recieves: 2 ETE subtrees and 2 node names
    returns: an ETEtree with the 2 concatenated ETE subtrees
    '''
    t_temp = Tree()  # for concatenation of both subtrees ahead, to avoid polytomy
    t_temp.add_child(prune_node_cp)
    t_curr = t_cp_p.copy()
    rgft_node_cp = t_curr & rgft_name  # locate the node in the copied subtree
    new_branch_length = rgft_node_cp.dist / 2
    
    rgft_loc = rgft_node_cp.up
    rgft_node_cp.detach()
    t_temp.add_child(rgft_node_cp, dist=new_branch_length)
    t_temp.name = nname
    rgft_loc.add_child(t_temp, dist=new_branch_length)  # regrafting
    if nname == "ROOT_LIKE":  # (4)
        t_temp.delete()
        # preserve the name of the root node, as this is a REAL node in this case
        preserve = True
    return t_curr, preserve


def SPR_generator(t):
    '''
    (1) avoid moves to a branch in the pruned subtree
    (2) avoid moves to the branch of a sibiling (retaining the same topology)
    (3) avoid moves to the branch leading to the parent of the (retaining the same topology)
    (4) handle the automatic "ROOTLIKE" node of ETE trees -
        if the PRUNE location is around the rootlike --> delete the "ROOT_LIKE" subnode when pasting in dest,
        and preserve the name of the REAL rootnode when converting back to newick
    '''
    
    prune_name = "N12"
    rgft_name = "N9"
    all_nodes_names = [n.name for n in t.get_descendants()]
    
    # perform pruning
    nname, subtree1, subtree2 = prune_branch(t,
                                             prune_name)  # subtree1 is the pruned subtree. subtree2 is the remaining subtree
    if nname == (t & rgft_name).up.name:  # captures (2)
        pass
    # continue
    remaining_nodes_names = [n.name for n in subtree2.get_descendants()]
    if rgft_name in remaining_nodes_names:  # captures (1) and (3)
        neighbor_tree, preserve = regraft_branch(subtree2, subtree1, rgft_name, nname)
        if preserve:  # namely, (4) is True
            neighbor_tree_str = neighbor_tree.write(format=1, format_root_node=True)
        else:
            neighbor_tree_str = neighbor_tree.write(format=1)
        print(neighbor_tree_str)
        print(Tree(newick=neighbor_tree_str, format=1).get_ascii(show_internal=True))
    else:
        print("Forbidden move: {} --> {} for this tree !".format(prune_name, rgft_name))
    # continue
    
    return




from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.Consensus import *
from Bio.Phylo.TreeConstruction import *
from ete3 import Tree

def generate_bootstrap_trees(dirpath, algo, nbootrees):
    '''
    algo could be either ‘nj’ or ‘upgma’
    '''
    msa = AlignIO.read(dirpath + "masked_species_real_msa.phy", format="phylip")
    calculator = DistanceCalculator('identity')
    constructor = DistanceTreeConstructor(calculator, algo)
    trees = bootstrap_trees(msa, nbootrees, constructor)
    
    return trees


def extract_branch_score(support_tree, branch_name):
    for nb in support_tree.find_clades():   # todo: replace with a function that locates a specific branch directly!
        nname = nb.name
        if not nname:
            continue
            
        if nname == branch_name:
            bsupport = nb.confidence
            bstrap_val = 0 if not bsupport else bsupport / 100
        
    return bstrap_val

    
from io import StringIO
def map_bootstraped(tree, trees):
    '''
    :param tree: a newick string
    :param trees: the bootstrap trees object returned by Biopython function
    :param pname: the name of the pruned_branch in the STARTING tree
    :param rname: the name of the rgrft_branch in the STARTING tree
    :return: the bootstrap value of both the pruned and the regraft branches
    '''
    
    support_tree = get_support(Phylo.read(StringIO(tree), "newick"), trees, NBOOTREES)
    print(support_tree)
    
    #Phylo.write(support_tree, dirpath + "test_support_tree_{}_biopython.txt".format((algo)), "newick")
    # t = Tree(dirpath + "test_support_tree_{}_biopython.txt".format((algo)), format=1)
    # print(t.get_ascii(show_internal=True))
    
    return support_tree




if __name__ == '__main__':
    '''
    tree_str_with_internal_names = "(((((Sp002:2.9e-07,Sp006:0.00173012)N8:0.00248844,Sp005:0.00978979)N6:0.0136241,Sp004:0.0388109)N4:0.0316201,Sp000:0.045075)N1:0.0837587,Sp003:0.00308638,Sp001:0.0180697);"
    tree_str_with_internal_names = "((0011:0.1,0012:0.3)N1:0.0625,0027:0.1,(0017:0.1,(0018:0.2,(0029:0.1,((0008:0.1,0006:0.1)N12:0.025,(0002:0.05,0031:0.1)N13:0.025)N11:0.0125)N9:0.15)N7:0.0125)N3:0.03125);"
    tree_str_with_internal_names = "(0008:0.0375,0006:0.1,((0031:0.05,0002:0.1)N4:0.025,(0029:0.1,(0018:0.05,((0027:0.1,(0011:0.1,0012:0.2)N15:0.1)N12:0.05,0017:0.3)N11:0.025)N9:0.025)N5:0.2)N3:0.01875);"
    t = Tree(newick=tree_str_with_internal_names, format=1)
    t.get_tree_root().name = "ROOT_LIKE"
    print(t.get_ascii(show_internal=True))
    SPR_generator(t)
    '''

    NBOOTREES = 300
    ALGO = 'nj'
    # todo: (1) check what are the internal node names in the resulting tree, namely how can I know what to lookup based on the prune and rgft names --> I need the prune.up.name & the orig rgft_name
    # todo: (2) separate the names i lookup in the starting tree (the original cut&paste names), and in the resulting tree (based on todo1 findings)
    
    df = pd.read_csv("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/summary_files/learning_all_moves_step1_test_new_features.csv", nrows=200)
    grouped_df_by_ds = df.groupby("path", sort=False)
    print(len(grouped_df_by_ds))
    for group_id, df_by_ds in grouped_df_by_ds:
        dirpath = df_by_ds["path"].values[0]
        print(dirpath)
        with open(dirpath + 'masked_species_real_msa.phy_phyml_tree_bionj.txt', 'r') as tree_fpr:
            startingtree = tree_fpr.read().strip()
            
        trees = generate_bootstrap_trees(dirpath, algo=ALGO, nbootrees=NBOOTREES)
        support_startingtree = map_bootstraped(startingtree, trees)
        
        dfr = pd.read_csv(dirpath + "newicks_step1.csv")
        for i, row in dfr.iterrows():
            for treetype in ['starting', 'resulting']:
                if "subtree" in row["rgft_name"]:
                    continue
                
                pname = row["prune_name"]
                if treetype == 'starting':
                    rname = row["rgft_name"]
                else:
                    tree = row["newick"]
                    prune_node = (Tree(tree, format=1)) & rname
                    if prune_node.up:  # todo: replace with the respective function (of .up) in biopython
                        rname = (prune_node.up).name
                        if not rname:
                            rname = "ROOT_LIKE"
                    else:
                        print("**************** 'up' does not exist for this node ***********")
                if "Sp" in pname and "Sp" in rname:
                    bstrap_prune, bstrap_rgft = 100, 100
                else:
                    support_tree = support_startingtree if treetype == 'starting' else map_bootstraped(tree, trees)  # if treetype=='resulting'
                    bstrap_prune = 100 if "Sp" in pname else extract_branch_score(support_tree, pname)
                    bstrap_rgft = 100 if "Sp" in rname else extract_branch_score(support_tree, rname)
                
                print(treetype)
                print(pname + ":", bstrap_prune)
                print(rname + ":", bstrap_rgft)
                print("########")

                df.loc[(df["path"] == dirpath) & (df["prune_name"] == pname) & (df["rgft_name"] == rname), "bstrap_{}_prune_{}".format(ALGO, treetype)] = bstrap_prune
                df.loc[(df["path"] == dirpath) & (df["prune_name"] == pname) & (df["rgft_name"] == rname), "bstrap_{}_rgft_{}".format(ALGO, treetype)] = bstrap_rgft
                
    df = df.head(200)
    print(len(df))
    print(df)
    df.to_csv("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/summary_files/learning_all_moves_step1_test_new_features_test.csv")
    df = df.dropna()
    print(len(df))
    
