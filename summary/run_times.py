import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")

# from defs import *
from utils.tree_functions import calc_leaves_features



def collect_features_runtime(ds_path, step_number):
    from pandas import read_csv
    from time import time
    TREES_PER_DS = "{}newicks_step{}.csv"
    PHYML_TREE_FILENAME = "masked_species_real_msa.phy" + "_phyml_tree_{0}.txt"

    dfr = read_csv(TREES_PER_DS.format(ds_path, step_number), index_col=0)
    '''
	from ete3 import Tree
	##print(Tree(ds_path + PHYML_TREE_FILENAME.format("bionj"), format=1).get_ascii(show_internal=True))
	from ete3 import TreeStyle
	ts = TreeStyle()
	ts.show_leaf_name = True
	ts.show_branch_length = True

	t = Tree(ds_path + PHYML_TREE_FILENAME.format("bionj"), format=1)
	##t.render(DATA_PATH+ "blah.png",tree_style=ts)
	#print(t.get_ascii(show_internal=True))
	#print("###############")
	#'''
    start_time = time()
    calc_leaves_features(ds_path + PHYML_TREE_FILENAME.format("bionj"), "prune")
    for i, row in dfr.iterrows():
        tree = row["newick"]
        if row["rgft_name"] == "subtree2":  # namely the remaining subtree
            calc_leaves_features(tree, "rgft")  # msa will be truncated within the function
        #if not "subtree" in row["rgft_name"]:
        #    start_time = time()
        #    calc_leaves_features(tree, "res", rgft_node_name=row["rgft_name"])

    res = time() - start_time
    return res


def extract_runtime(statsfile):
    try:
        with open(statsfile) as fpr:
            content = fpr.read()
        # runtime = re.search("Time used:\s+(\d+)h(\d+)m(.*?)s", content)
        # hours, minutes, seconds = runtime.group(1), runtime.group(2), runtime.group(3)
        # print(hours,minutes,seconds)
        runtime = re.search("Time used:\s+(\d+)h(\d+)m(.*?)s\s+\((.*)\s+seconds\)", content).group(4)
    except:
        print("Error with:", statsfile)
        return

    return int(runtime)


def calc_ll_runtime(dirname):
    secs = 0
    for prune_dir in os.listdir(dirname):
        prune_path = dirname + prune_dir
        for rgft_dir in os.listdir(prune_path):
            if not "subtree" in rgft_dir:
                full_dirpath = SEP.join([prune_path, rgft_dir, ""])
                # msa = full_dirpath + MSA_PHYLIP_FILENAME
                # res = call_phyml(full_dirpath, REARRANGEMENTS_NAME, msa, True, -1, "br ",cpmsa=False)
                phyml_stats = PHYML_STATS_FILENAME.format('br')
                res_secs = extract_runtime(full_dirpath + phyml_stats)
                if res_secs:
                    secs += res_secs
                else:
                    print(full_dirpath)
    return secs


if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser(description='arrange data for learning and implement learning algo')
    parser.add_argument('--ll_or_features', '-t', required=True)  # ll OR features
    args = parser.parse_args()

    DATA_PATH = r"C:\Users\ItayMNB3\Dropbox\PhyloAI\data\training_datasets\\"
    # DATA_PATH = "/groups/itay_mayrose/danaazouri/PhyAI/DBset2/data/training_datasets/"
    # dirname = DATA_PATH + "example4/rearrangements/"
    dirname = DATA_PATH + "test/rearrangements/"

    if args.ll_or_features == "ll":
        res = calc_ll_runtime(dirname)
    if args.ll_or_features == "features":
        # res = collect_features_runtime(DATA_PATH + "example4/", "1")
        res = collect_features_runtime("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/data/training_datasets/exampleSphaero/", "1")

    print("--- %s seconds ---" % (str(res)))