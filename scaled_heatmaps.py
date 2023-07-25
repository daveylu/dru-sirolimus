import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.gridspec
import matplotlib.cm
import matplotlib.pyplot as plt

# metrics = ["TAT 4 wk", "TAT 20 wk", "Cmax", "Cmax4"]          # possible metrics to use
# markers = ["4", "20"]                                         # week of metrics taken

plt.rcParams['text.usetex'] = True                              # enables LaTeX


# metric to sort by: name of column in excel file to sort by
sort_via = "TAT 4 wk"

# flag: primary only or primary and secondary
primary_only = True

# flag: 4 week biomarkers or 20 week biomarkers
separate_biomarkers = True
markers = "20"

# flag: only view significant biomarkers
significant_biomarkers_only = False
significance_level = 0.1

# read in our datasets
val_df = pd.read_csv("modified_data/scaled_fold_change_no_outliers.csv")    # contains fold change data per patient and biomarker
pval_df = pd.read_csv("modified_data/rho_pvals.csv")            # contains p-values of biomarkers per metric

if(primary_only):
    val_df = val_df[val_df["in_primary"] == 1]                  # only use patients in the primary group by filtering out non-primary patients
    title_str = "Primary Group Only"
else:
    title_str = "Primary & Secondary Groups"


# sort and extract relevant p-value data
pval_df = pval_df[pval_df["exposure"] == sort_via]              # find p-vals corresponding to the metric we want to sort by

if(significant_biomarkers_only):                                # only want to look at biomarkers that are significant
    pval_df = pval_df[pval_df["pval"] < significance_level]     # keep rows with p-values below the level we choose
    biomarkers_to_keep = list(pval_df["biomarker"])             # get list of the biomarkers that match this criteria
    patient_info_df = val_df.iloc[:, 0:6]
    val_df = val_df.loc[:, biomarkers_to_keep]                  # only keep the data from the biomarkers that matched criteria
    val_df = pd.concat([patient_info_df, val_df], axis = 1)

pval_df = pval_df.sort_values(by = ["pval"]) # type: ignore     # sort p-vals from lowest to highest
p_vals = pval_df["pval"]                                        # extract pd.Series of sorted p-values
rhos = pval_df["rho"]                                           # extract pd.Series of rho values

# keep biomarkers as determined by the flag if we want to separate them: 4 week or 20 week
if(separate_biomarkers == True):
    patient_info_df = val_df.iloc[:, 0:6]                       # splitting up the dataframe into the patient info and the actual data
    val_df = val_df.iloc[:, 6:]
    markers = "." + markers                                     # minor string modification
    cols_to_keep = val_df.columns.str.endswith(markers)         # creates boolean array corresponding to which columns are from the correct week
    val_df = val_df.loc[:, cols_to_keep]                        # grabs dataframe with the correct columns only
    val_df = pd.concat([patient_info_df, val_df], axis = 1)

# sort and extract relevant biomarker and fold change data
num_patients, num_cols = val_df.shape                           # find number of patients (rows) and number of columns
num_biomarkers = num_cols - 6                                   # num_biomarkers = num_cols - number of non-biomarker columns
ids = val_df["Patient ID"]                                      # extract pd.Series of patient IDs (Patient ID column)
col_to_sort = list(val_df.columns).index(sort_via)              # find column that we want to use to sort by
indices = [0, col_to_sort] + list(range(6, num_cols))           # get indices of columns to take the patient IDs, metric to sort by, and all data points
val_df = val_df.iloc[:, indices]                                # grab those columns in our data, toss out the rest
val_df = val_df.fillna(0)                                       # replace NaNs (Not a Number) with 0 since sns.clustermap cannot have NaNs

val_df = val_df.sort_values(by = [sort_via])                    # sort by sort_via metric: lowest -> highest
exposure = val_df[sort_via]                                     # extract pd.Series of exposure metrics

df_plot = val_df.transpose()                                    # transpose our dataframe: a column is now a patient
df_plot.columns = df_plot.iloc[0]                               # set our column labels to be our patient IDs
df_plot = df_plot.iloc[2:]                                      # get rid of the patient IDs row and exposure metric rows

# create objects for coloring patient exposure, p-values, and rho values
exposure_cmap = sns.color_palette("flare", as_cmap = True)                                                      # create exposure colormap (cmap) object
exposure_norm = matplotlib.colors.Normalize(vmin = exposure.min(), vmax = exposure.max())                       # create exposure normalizer 
exposure_mapper = matplotlib.cm.ScalarMappable(norm = exposure_norm, cmap = exposure_cmap) #type: ignore        # create mapper which takes in a patient's exposure and outputs a color
exposure_colors = [exposure_mapper.to_rgba(pat) for pat in exposure]                                            # create list of colors, from lowest exposure to highest

# create color palette series: df_plot.columns contains sorted patient IDs, used as dict keys to index into the palette
patient_color_palette = pd.Series(data = exposure_colors, index = df_plot.columns, name = f"{sort_via} Exposure")

pval_cmap = sns.color_palette("viridis", as_cmap = True)
pval_norm = matplotlib.colors.Normalize(vmin = 0, vmax = 1)
pval_mapper = matplotlib.cm.ScalarMappable(norm = pval_norm, cmap = pval_cmap)  #type: ignore
pval_colors = [pval_mapper.to_rgba(p_val) for p_val in p_vals]

rho_cmap = sns.color_palette("BrBG", as_cmap = True)
rho_norm = matplotlib.colors.Normalize(vmin = -1, vmax = 1)
rho_mapper = matplotlib.cm.ScalarMappable(norm = rho_norm, cmap = rho_cmap)     #type: ignore
rho_colors = [rho_mapper.to_rgba(rho) for rho in rhos]

# create a pd.DataFrame that contains the corresponding p-values and rho values for each biomarker
row_color_palette = pd.DataFrame({"p-values": pval_colors,
                                 r"Spearman's $\rho$": rho_colors},
                                index = pval_df["biomarker"])

# plots our clustered heatmap, using our dataframe from above and clustering via rows ONLY
# patients have increasing exposure from left to right: lowest -> highest exposure according to set metric
# clustering: hierarchal clustering w/ euclidean distance used (default)
# NOTE: scipy must be installed for sns.clustermap to work
clustermap = sns.clustermap(df_plot, metric = "seuclidean", row_cluster = True, col_cluster = False, cmap = "mako",
                            row_colors = row_color_palette, col_colors = patient_color_palette,
                             linewidth = 0.5, xticklabels = True, yticklabels = True,
                            cbar_kws = {"label": "Scaled Fold Change"}, cbar_pos = (0.89, 0.2, 0.015, 0.6))

# plots colorbars for p-values and biomarkers onto the figure
exposure_ax = clustermap.figure.add_axes([0.2, 0.9, 0.6, 0.025])                                # add a new axes onto the figure for the colorbar
exposure_cbar = plt.colorbar(mappable = exposure_mapper, cax = exposure_ax,                     # create a new color bar using the mapper
                             orientation = "horizontal", label = f"{sort_via} Exposure")        # corresponding to that colorbar

pval_ax = clustermap.figure.add_axes([0.94, 0.525, 0.015, 0.275])
pval_cbar = plt.colorbar(mappable = pval_mapper, cax = pval_ax,
                         orientation = "vertical", label = f"{sort_via} p-values")

rho_ax = clustermap.figure.add_axes([0.94, 0.2, 0.015, 0.275])
rho_cbar = plt.colorbar(mappable = rho_mapper, cax = rho_ax, orientation = "vertical",
                        label = rf"{sort_via} Spearman's $\rho$ values")

# add title to the figure
clustermap.figure.suptitle(f"Scaled Biomarker Fold Change\n{title_str}", fontsize = "xx-large")

plt.show()  # displays the plot