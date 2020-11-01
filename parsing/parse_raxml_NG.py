import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")

from defs import *
from utils import msa_functions



def parse_raxmlNG_output(res_filepath):

    try:
        with open(res_filepath) as fpr:
            content = fpr.read()
        res_dict = parse_raxmlNG_content(content)
    except:
        print("Error with:", res_filepath)
        return

    return res_dict


def parse_raxmlNG_content(content):
    """
	:return: dictionary with the attributes - string typed. if parameter was not estimated, empty string
	"""
    res_dict = dict.fromkeys(["ll", "pInv", "gamma",
                              "fA", "fC", "fG", "fT",
                              "subAC", "subAG", "subAT", "subCG", "subCT", "subGT",
                              "time"], "")

    # likelihood
    ll_re = re.search("Final LogLikelihood:\s+(.*)", content)
    if ll_re:
        res_dict["ll"] = ll_re.group(1).strip()
    elif re.search("BL opt converged to a worse likelihood score by", content) or re.search("failed", content):
        ll_ini = re.search("initial LogLikelihood:\s+(.*)", content)
        if ll_ini:
            res_dict["ll"] = ll_ini.group(1).strip()
    else:
        res_dict["ll"] = 'unknown raxml-ng error, check "parse_raxmlNG_content" function'


    # gamma (alpha parameter) and proportion of invariant sites
    gamma_regex = re.search("alpha:\s+(\d+\.?\d*)\s+", content)
    pinv_regex = re.search("P-inv.*:\s+(\d+\.?\d*)", content)
    if gamma_regex:
        res_dict['gamma'] = gamma_regex.group(1).strip()
    if pinv_regex:
        res_dict['pInv'] = pinv_regex.group(1).strip()

    # Nucleotides frequencies
    nucs_freq = re.search("Base frequencies.*?:\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)", content)
    if nucs_freq:
        for i,nuc in enumerate("ACGT"):
            res_dict["f" + nuc] = nucs_freq.group(i+1).strip()

    # substitution frequencies
    subs_freq = re.search("Substitution rates.*:\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)", content)
    if subs_freq:
        for i,nuc_pair in enumerate(["AC", "AG", "AT", "CG", "CT", "GT"]):  # todo: make sure order
            res_dict["sub" + nuc_pair] = subs_freq.group(i+1).strip()

    # Elapsed time of raxml-ng optimization
    rtime = re.search("Elapsed time:\s+(\d+\.?\d*)\s+seconds", content)
    if rtime:
        res_dict["time"] = rtime.group(1).strip()
    else:
        res_dict["time"] = 'no ll opt_no time'

    return res_dict



if __name__ == '__main__':
    pass
    # test here
    res_dict = parse_raxmlNG_output("/groups/itay_mayrose/danaazouri/PhyAI/DBset2/data_test/real_msa3.phy.raxml.log")
    print(res_dict)
    print(res_dict["time"])
