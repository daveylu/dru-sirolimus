import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats

# flag: primary only or primary and secondary
primary_only = True

#flag: 4 week biomarker data or 20 week biomarker data
separate_weeks = True
week = "20"

if(primary_only):
    title_str = "Primary Group"
else:
    title_str = "Primary and Secondary Group"

plt.rcParams['text.usetex'] = True                                                      # enables LaTeX

df = pd.read_csv("modified_data/fold_change_no_outliers.csv")                           # load dataset

if(primary_only):                                                                       # keep only patients in primary group if flag is set
    df = df[df["in_primary"] == 1]

biomarkers_df = df.iloc[:, 6:]                                                          # keep only the biomarker data

if(separate_weeks):
    cols_to_keep = biomarkers_df.columns.str.endswith(week)
    biomarkers_df = biomarkers_df.loc[:, cols_to_keep]

biomarkers = list(biomarkers_df.columns) #type: ignore
num_biomarkers = biomarkers_df.shape[1]  #type: ignore

corr_mat, p_vals = scipy.stats.spearmanr(biomarkers_df, nan_policy = "omit")            # calculate Spearman's rho          # type: ignore
corr_df = pd.DataFrame(data = corr_mat, index = biomarkers, columns = biomarkers)       # stick the resulting np.ndArray into a pd.DataFrame

# create a mask to block out the upper triangular portion of the correlation matrix
mask = np.ones((num_biomarkers, num_biomarkers))
mask = np.triu(mask)

fig = sns.clustermap(corr_df, row_cluster = True, col_cluster = True, cmap = "seismic",
                   xticklabels = True, yticklabels = True, linewidth = 0.5,
                   cbar_kws = {"label": r"Spearman's $\rho$"}, vmin = -1, vmax = 1,
                   cbar_pos = [0.925, 0.2, 0.025, 0.6], figsize = (13, 14))

# below section performs the masking AFTER the clustering is done
# taken from: https://stackoverflow.com/questions/67879908/lower-triangle-mask-with-seaborn-clustermap
values = fig.ax_heatmap.collections[0].get_array().reshape(corr_df.shape)
new_values = np.ma.array(values, mask=mask)
fig.ax_heatmap.collections[0].set_array(new_values)

fig.figure.text(0.1, 0.9, title_str, fontsize = "xx-large",                            # add a title to the plot
                horizontalalignment = "center", verticalalignment = "center")

plt.show()