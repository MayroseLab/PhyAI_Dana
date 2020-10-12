import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")
import warnings

from defs import *
from Bio import SeqIO
#import networkx, pylab

#PATH = "/groups/itay_mayrose/danaazouri/PhyAI/DBset2/data/training_datasets/exampleSphaero/real_msa"
#SeqIO.convert(PATH + '.nex', "nexus", PATH + '.phy', PHYLIP_FORMAT)

def draw_phytree():
    tree = "(((Sp004:0.00978979,Sp003:0.0388109)N2:0.0316201,Sp002:0.045075)N1:0.0837587,Sp001:0.00308638)ROOT;"
    t = Tree(tree, format=1)
    print(t.get_ascii(show_internal=True), "\n###########\n")
    '''
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
    for tn in [t]:
        for l in tn.iter_leaves():
            l.img_style = nst_general

    nst1 = NodeStyle()
    # nst1["bgcolor"] = "LightSteelBlue"n #"Moccasin"
    nst1["fgcolor"] = "red"
    nst1["size"] = 5

    nst2 = NodeStyle()
    nst2["fgcolor"] = "blue"
    nst2["size"] = 5
    '''
    t.render(r"C:\Users\ItayMNB3\Dropbox\PhyloAI\\" + 'trees.png') #, tree_style=ts)




if __name__ == '__main__':
    #draw_phytree()

    df = pd.read_csv(SUMMARY_FILES_DIR + "with_preds_merged_20_1_5850_ytransformed_exp.csv", nrows=1000)

    #corr_path = "/groups/itay_mayrose/danaazouri/PhyAI/DBset2/data/training_datasets/24518/"
    corr_path = "/groups/itay_mayrose/danaazouri/PhyAI/DBset2/data/training_datasets/38663/"
    df = df[df["path"] == corr_path]
    #print(df[["d_ll_merged", "pred"]].corr(method='spearman'))
    #print(df)

    #sns.heatmap(corr, annot=True)
    plt.scatter(df["d_ll_merged"].rank(), df["pred"].rank())
    plt.xlabel("Rank/naccording to the log-likelihood", labelpad=13)
    plt.ylabel("Rank/naccording to the model", labelpad=13)
    plt.tight_layout()
    plt.savefig(SUMMARY_FILES_DIR + "corr.png")
    #plt.show()