import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")

from defs import *
from utils import msa_functions


def parse_raxmlNG_output(res_filepath):
    """
	:return: dictionary with the attributes - string typed. if parameter was not estimated, empty string
	"""
    res_dict = dict.fromkeys(["ll", "pInv", "gamma",
                              "fA", "fC", "fG", "fT",
                              "subAC", "subAG", "subAT", "subCG", "subCT", "subGT"], "")

    try:
        with open(res_filepath) as fpr:
            content = fpr.read()

        # likelihood
        res_dict["ll"] = re.search("Final LogLikelihood:\s+(.*)", content).group(1).strip()
        if not res_dict['ll'] and re.search("BL opt converged to a worse", content):
            res_dict["ll"] = re.search("initial LogLikelihood:\s+(.*)", content).group(1).strip()

        # gamma (alpha parameter) and proportion of invariant sites
        gamma_regex = re.search("alpha:\s+(\d+\.?\d*)\s+", content)
        pinv_regex = re.search("P-inv.*:\s+(\d+\.?\d*)", content)
        if gamma_regex:
            res_dict['gamma'] = gamma_regex.group(1).strip()
        if pinv_regex:
            res_dict['pInv'] = pinv_regex.group(1).strip()

        # Nucleotides frequencies
        nucs_freq = re.search("Base frequencies.*?:\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)", content)
        for i,nuc in enumerate("ACGT"):
            res_dict["f" + nuc] = nucs_freq.group(i+1).strip()

        # substitution frequencies
        subs_freq = re.search("Substitution rates.*:\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)", content)
        for i,nuc_pair in enumerate(["AC", "AG", "AT", "CG", "CT", "GT"]):  # todo: make sure order
            res_dict["sub" + nuc_pair] = subs_freq.group(i+1).strip()

    except:
        print("Error with:", res_filepath)
        return

    return res_dict



if __name__ == '__main__':
    pass
    # test here
   #res_dict = parse_raxmlNG_output("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/data_test/real_msa4.phy.raxml.log")
   # print(res_dict)
