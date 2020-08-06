dirpath = r'C:\Users\ItayMNB3\Dropbox\courses\python_bio\2020\projects\Dana\submissions\\'

#####################

import pandas as pd
import numpy as np
import scipy
import os
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
import math
from scipy.spatial import distance_matrix
from scipy import stats
from scipy.stats import zscore


# Functions for script ##########################################################################
# Z score function
def feature_stn(table):
    norm_table = table.apply(zscore)
    return norm_table


# Regression functions
def reg(feature_vec, lipo_vec, cova):
    df = pd.concat([feature_vec, lipo_vec, cova["Age"], cova["BMI"]], axis=1)
    df.columns = ["feature_obs", "lipo_obs", "Age", "BMI"]
    df = df.dropna(axis=0)
    fit = LinearRegression().fit(df.loc[:, ["lipo_obs", "Age", "BMI"]], df["feature_obs"])
    return fit.coef_[0]


def lp_connectivity_func(features, lipos, cova, groups_dict):
    lp_dfs_dict = {}
    for group_name, group in groups_dict.items():
        group_features = features.loc[group, :]
        group_lipos = lipos.loc[group, :]
        group_cova = cova.loc[group, :]
        lp_table = group_features.apply(
            lambda fea_col: group_lipos.apply(lambda lipo_col: reg(fea_col, lipo_col, group_cova), axis=0), axis=0)
        lp_dfs_dict[group_name] = lp_table
    return lp_dfs_dict


# Rotate x,y around xo,yo by theta (rad)
def rotate(x, y, xo, yo, theta):
    theta = theta * math.pi / 180
    cen_x = (x.max() + x.min()) / 2
    cen_y = (y.max() + y.min()) / 2
    x = np.subtract(x, cen_x)
    y = np.subtract(y, cen_y)
    xr = -1 * (np.cos(theta) * np.subtract(x, xo) - np.sin(theta) * np.subtract(y, yo) + xo)
    yr = np.sin(theta) * np.subtract(x, xo) + np.cos(theta) * np.subtract(y, yo) + yo
    return pd.DataFrame({"PC1": xr, "PC2": yr})


# Calculate tne angel of the feature in the cycle
def calc_angel(pca_table):
    pca_table.columns = ["PC1", "PC2"]
    angel = np.arctan2(pca_table["PC2"], pca_table["PC1"]) * (180 / math.pi)
    angel = angel.apply(lambda x: x + 360 if x < 0 else x)
    return angel


# Plot PCA graph for each group, all group are colored by the healthy females group
# Save the plot in the path folder
def pca_for_groups(lp_dfs_dict, females_pca_model):
    pcas_dict = {}
    for i, group in enumerate(lp_dfs_dict.keys()):
        pca = pd.DataFrame(females_pca_model.transform(feature_stn(lp_dfs_dict[group]).T))
        pca.columns = ["PC1", "PC2"]
        pca.index = lp_dfs_dict[group].T.index
        pca = rotate(pca["PC1"], pca["PC2"], 0, 0, 110)
        pca["angel"] = calc_angel(pca)
        pcas_dict[group] = pca
    fig, axs = plt.subplots(1, 4, figsize=(28, 7))
    for i, group in enumerate(pcas_dict.keys()):
        axs[i].scatter(x="PC1", y="PC2", c=pcas_dict["fe_no_MetS_no_AT"]["angel"], cmap="hsv", data=pcas_dict[group])
        axs[i].title.set_text(group)
        axs[i].set_xlabel("PC1")
        axs[i].set_ylabel("PC2")
        axs[i].axhline(y=0, color="black")
        axs[i].axvline(x=0, color="black")
    plt.savefig(os.path.join(folder_path, "PCA_graphs.PNG"))
    plt.show()


# Plot lp connectivity cluster maps for each group
# Save the plot in the path folder
def plot_cluster_maps(lp_dict):
    for group_name, group in lp_dict.items():
        g = sns.clustermap(lp_dict[group_name], figsize=(35, 18), cmap="bwr")
        g.fig.suptitle(group_name + ' cluster map', fontsize="xx-large", y=0.9)
        ax = g.ax_heatmap
        ax.set_xlabel("Clinical features", fontsize="xx-large")
        ax.set_ylabel("Lipoproteins", fontsize="xx-large")
        plt.savefig(os.path.join(folder_path, group_name + "_cluster_map.PNG"))
        plt.show()


# Use the function for analysis ##########################################################################

# Change this path the load tables
folder_path = dirpath

# read tables
lipoNorm = pd.read_csv(os.path.join(folder_path, "lipoNorm_fem.csv"), index_col=0)
covaNorm = pd.read_csv(os.path.join(folder_path, "covaNorm_fem.csv"), index_col=0)
factors = pd.read_csv(os.path.join(folder_path, "factors_fem.csv"), index_col=0)
factorsNorm = pd.read_csv(os.path.join(folder_path, "FactorsNorm_fem.csv"), index_col=0)
subjects_anott = pd.read_csv(os.path.join(folder_path, "subjects_anott_fem.csv"), index_col=0)

subjects_anott["atherosclerosis"] = factors["Plaques aanwezig"]
factorsNorm = factorsNorm.drop(["PWV_direct_def", "Numberofplaques"], axis=1)

# AT = Atherosclerosis (Arterial stenosis)
# MetS = Metabolic syndrome

# Create a Index vector for each group subjects:
# Females healthy,
# Females without MetS with AT
# Females with MetS without AT
# Females with MetS with AT
fem_groups_dict = {"fe_no_MetS_no_AT": subjects_anott[(subjects_anott["MetS-IDFcriteria"] == 0) &
                                                      (subjects_anott["atherosclerosis"] == 0)].index,
                   "fe_no_MetS_AT": subjects_anott[(subjects_anott["MetS-IDFcriteria"] == 0) &
                                                   (subjects_anott["atherosclerosis"] == 1)].index,
                   "fe_MetS_no_AT": subjects_anott[(subjects_anott["MetS-IDFcriteria"] == 1) &
                                                   (subjects_anott["atherosclerosis"] == 0)].index,
                   "fe_MetS_AT": subjects_anott[(subjects_anott["MetS-IDFcriteria"] == 1) &
                                                (subjects_anott["atherosclerosis"] == 1)].index}


# Use the regression function the calculate the lp connectivity of the clinical features in each group
lp_connect_dict = lp_connectivity_func(factorsNorm, lipoNorm, covaNorm, fem_groups_dict)

# Use the regression function the calculate the lp connectivity of the lipoproteins of the healthy group
lp_of_lps_healthy = lp_connectivity_func(lipoNorm, lipoNorm, covaNorm,
                                          {"fe_no_MetS_no_AT": fem_groups_dict["fe_no_MetS_no_AT"]})

# Calculate a PCA model of the healthy subjects
lps_pca_model_fem = PCA(n_components=2).fit(feature_stn(lp_of_lps_healthy["fe_no_MetS_no_AT"]).T)

pca_lipo_fem = pd.DataFrame(lps_pca_model_fem.transform(feature_stn(lp_of_lps_healthy["fe_no_MetS_no_AT"]).T))
pca_lipo_fem.columns = ["PC1", "PC2"]
pca_lipo_fem.index = lipoNorm.T.index

pca_lipo_fem = rotate(pca_lipo_fem["PC1"], pca_lipo_fem["PC2"], 0, 0, 110)
pca_lipo_fem["angel"] = calc_angel(pca_lipo_fem)

# Use the plotting function the plot PCA for each group
pca_for_groups(lp_connect_dict, lps_pca_model_fem)

# Use the clustermap function the plot clustermap for each group
plot_cluster_maps(lp_connect_dict)
