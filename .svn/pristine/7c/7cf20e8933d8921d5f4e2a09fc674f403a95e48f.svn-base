import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")

from defs import *
from utils import msa_functions


def parse_phyml_stats_output(msa_filepath, stats_filepath):
    """
	:return: dictionary with the attributes - string typed. if parameter was not estimated, empty string
	"""
    res_dict = dict.fromkeys(["ntaxa", "nchars", "ll",
                              "fA", "fC", "fG", "fT",
                              "subAC", "subAG", "subAT", "subCG", "subCT", "subGT",
                              "pInv", "gamma",
                              "path"], "")
    
    if msa_filepath:
        res_dict['ntaxa'], res_dict['nchars'] = (str(x) for x in msa_functions.get_msa_properties(msa_functions.get_msa_from_file(msa_filepath)))
    
    res_dict["path"] = stats_filepath
    try:
        with open(stats_filepath) as fpr:
            content = fpr.read()
        
        # likelihood
        res_dict["ll"] = re.search("Log-likelihood:\s+(.*)", content).group(1).strip()
        
        # gamma (alpha parameter) and proportion of invariant sites
        gamma_regex = re.search("Gamma shape parameter:\s+(.*)", content)
        pinv_regex = re.search("Proportion of invariant:\s+(.*)", content)
        if gamma_regex:
            res_dict['gamma'] = gamma_regex.group(1).strip()
        if pinv_regex:
            res_dict['pInv'] = pinv_regex.group(1).strip()
        
        # Nucleotides frequencies
        for nuc in "ACGT":
            nuc_freq = re.search("  - f\(" + nuc + "\)\= (.*)", content).group(1).strip()
            res_dict["f" + nuc] = nuc_freq
        
        # substitution frequencies
        for nuc1 in "ACGT":
            for nuc2 in "ACGT":
                if nuc1 < nuc2:
                    nuc_freq = re.search(nuc1 + " <-> " + nuc2 + "(.*)", content).group(1).strip()
                    res_dict["sub" + nuc1 + nuc2] = nuc_freq
    except:
        print("Error with:", res_dict["path"], res_dict["ntaxa"], res_dict["nchars"])
        return
    return res_dict


if __name__ == '__main__':
    pass
    # test here
