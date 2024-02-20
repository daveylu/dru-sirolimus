import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.gridspec
import matplotlib.cm
import matplotlib.pyplot as plt

# metrics = ["TAT 4 wk", "TAT 20 wk", "Cmax", "Cmax4"]              # possible metrics to use
# markers = ["4", "20"]                                             # week of metrics taken

plt.rcParams['text.usetex'] = True                                  # enables LaTeX
fold_change_path = "modified_data/minmax_scaled_fold_change_no_outliers.csv"

def display_clustermap(sort_via: str, primary_only: bool = False, separate_biomarkers: bool = False, markers: str | None = None,
                       significant_biomarkers_only: bool = False, max_pval: float | None = None) -> None:
    """
    Displays a clustered heatmap from the dataset
    
    Parameters
    ----------
    sort_via: str
        metric to sort by
    primary_only: bool
        plot data from primary group only or primary and secondary groups
    separate_biomarkers: bool
        plot data from only one week of biomarkers or both weeks of biomarkers
    markers: str | None
        choose which week of biomarkers to plot if separate_biomarkers == True
        possible values: "4" or "20"
    significant_biomarkers_only: bool
        choose whether to only view biomarkers with a p-value below a desired value or not
    max_pval: float | None
        choose what is the maximum p-value to display if significant_biomarkers_only == True
        0 <= max_pval <= 1
    """

    # read in our datasets
    val_df = pd.read_csv(fold_change_path)                              # contains fold change data per patient and biomarker
    pval_df = pd.read_csv("modified_data/rho_pvals.csv")                # contains p-values of biomarkers per metric

    if(primary_only):
        val_df = val_df[val_df["in_primary"] == 1]                      # only use patients in the primary group by filtering out non-primary patients
        title_str = "Primary Group Only"
    else:
        title_str = "Primary and Secondary Groups"


    # sort and extract relevant p-value data
    pval_df = pval_df[pval_df["exposure"] == sort_via]                  # find p-vals corresponding to the metric we want to sort by

    if(significant_biomarkers_only):                                    # only want to look at biomarkers that are significant
        pval_df = pval_df[pval_df["pval"] < max_pval]                   # keep rows with p-values below the level we choose
        biomarkers_to_keep = list(pval_df["biomarker"])                 # get list of the biomarkers that match this criteria
        patient_info_df = val_df.iloc[:, 0:6]
        val_df = val_df.loc[:, biomarkers_to_keep]                      # only keep the data from the biomarkers that matched criteria
        val_df = pd.concat([patient_info_df, val_df], axis = 1)

    pval_df = pval_df.sort_values(by = ["pval"])                        # sort p-vals from lowest to highest
    p_vals = pval_df["pval"]                                            # extract pd.Series of sorted p-values
    rhos = pval_df["rho"]                                               # extract pd.Series of rho values

    # keep biomarkers as determined by the flag if we want to separate them: 4 week or 20 week
    if(separate_biomarkers == True):
        patient_info_df = val_df.iloc[:, 0:6]                           # splitting up the dataframe into the patient info and the actual data
        val_df = val_df.iloc[:, 6:]
        markers = "." + markers                                         # minor string modification         #type: ignore
        cols_to_keep = val_df.columns.str.endswith(markers)             # creates boolean array corresponding to which columns are from the correct week    # type: ignore
        val_df = val_df.loc[:, cols_to_keep]                            # grabs dataframe with the correct columns only
        val_df = pd.concat([patient_info_df, val_df], axis = 1)

    # sort and extract relevant biomarker and fold change data
    num_patients, num_cols = val_df.shape                               # find number of patients (rows) and number of columns
    num_biomarkers = num_cols - 6                                       # num_biomarkers = num_cols - number of non-biomarker columns
    ids = val_df["Patient ID"]                                          # extract pd.Series of patient IDs (Patient ID column)
    col_to_sort = list(val_df.columns).index(sort_via)                  # find column that we want to use to sort by
    indices = [0, col_to_sort] + list(range(6, num_cols))               # get indices of columns to take the patient IDs, metric to sort by, and all data points
    val_df = val_df.iloc[:, indices]                                    # grab those columns in our data, toss out the rest
    val_df = val_df.fillna(0)                                           # replace NaNs (Not a Number) with 0 since sns.clustermap cannot have NaNs

    val_df = val_df.sort_values(by = [sort_via])                        # sort by sort_via metric: lowest -> highest
    exposure = val_df[sort_via]                                         # extract pd.Series of exposure metrics

    df_plot = val_df.transpose()                                        # transpose our dataframe: a column is now a patient
    df_plot.columns = df_plot.iloc[0]                                   # set our column labels to be our patient IDs
    df_plot = df_plot.iloc[2:]                                          # get rid of the patient IDs row and exposure metric rows

    # create objects for coloring patient exposure, p-values, and rho values
    exposure_cmap = sns.color_palette("flare", as_cmap = True)                                                      # create exposure colormap (cmap) object
    exposure_norm = matplotlib.colors.Normalize(vmin = exposure.min(), vmax = exposure.max())                       # create exposure normalizer 
    exposure_mapper = matplotlib.cm.ScalarMappable(norm = exposure_norm, cmap = exposure_cmap) #type: ignore        # create mapper which takes in a patient's exposure and outputs a color
    exposure_colors = [exposure_mapper.to_rgba(pat) for pat in exposure]                                            # create list of colors, from lowest exposure to highest

    # create color palette series: df_plot.columns contains sorted patient IDs, used as dict keys to index into the palette
    patient_color_palette = pd.Series(data = exposure_colors, index = df_plot.columns, name = f"{sort_via} Exposure")

    # create colors for p-values
    pval_cmap = sns.color_palette("viridis", as_cmap = True)
    if(significant_biomarkers_only):                                                                                # set vmax to maximum value allowed for p-values
        pval_norm = matplotlib.colors.Normalize(vmin = 0, vmax = max_pval)                                          # when looking at significant biomarkers only
    else:
        pval_norm = matplotlib.colors.Normalize(vmin = 0, vmax = 1)
    pval_mapper = matplotlib.cm.ScalarMappable(norm = pval_norm, cmap = pval_cmap)  #type: ignore
    pval_colors = [pval_mapper.to_rgba(p_val) for p_val in p_vals]

    # repeat above except for rho values
    rho_cmap = sns.color_palette("BrBG", as_cmap = True)
    rho_norm = matplotlib.colors.Normalize(vmin = -0.75, vmax = 0.75)
    rho_mapper = matplotlib.cm.ScalarMappable(norm = rho_norm, cmap = rho_cmap)     #type: ignore
    rho_colors = [rho_mapper.to_rgba(rho) for rho in rhos]

    # create a pd.DataFrame that contains the corresponding p-values and rho values for each biomarker
    row_color_palette = pd.DataFrame({"p-values": pval_colors,
                                    r"Spearman's $\rho$": rho_colors},
                                    index = pval_df["biomarker"])

    # plots our clustered heatmap, using our dataframe from above and clustering via rows ONLY
    # patients have increasing exposure from left to right: lowest -> highest exposure according to set metric
    # clustering: hierarchal clustering w/ standardized euclidean distance used
    # NOTE: scipy must be installed for sns.clustermap to work

    print(df_plot.to_numpy())

    clustermap = sns.clustermap(df_plot, metric = "seuclidean", row_cluster = True, col_cluster = False, linewidth = 0.5,
                                row_colors = row_color_palette, col_colors = patient_color_palette, cmap = "plasma",
                                xticklabels = True, yticklabels = True, vmin = 0, vmax = 1, center = 0.5, figsize=(20, 12),
                                cbar_kws = {"label": "Fold Change", "orientation": "horizontal"}, cbar_pos = (0.6, 0.9, 0.3, 0.025))

    if(not significant_biomarkers_only):        # reduce font size of biomarkers so they don't overlap
        clustermap.ax_heatmap.set_yticklabels(clustermap.ax_heatmap.get_ymajorticklabels(), fontsize = 6)

    # plots colorbars for p-values and biomarkers onto the figure
    exposure_ax = clustermap.figure.add_axes([0.2, 0.9, 0.3, 0.025])                                # add a new axes onto the figure for the colorbar
    exposure_cbar = plt.colorbar(mappable = exposure_mapper, cax = exposure_ax,                     # create a new color bar using the mapper
                                 orientation = "horizontal", label = f"{sort_via} Exposure")        # corresponding to that colorbar

    pval_ax = clustermap.figure.add_axes([0.05, 0.85, 0.02, 0.1])
    pval_cbar = plt.colorbar(mappable = pval_mapper, cax = pval_ax,
                            orientation = "vertical", label = f"{sort_via} p-values")

    rho_ax = clustermap.figure.add_axes([0.1, 0.85, 0.02, 0.1])
    rho_cbar = plt.colorbar(mappable = rho_mapper, cax = rho_ax, orientation = "vertical",
                            label = rf"{sort_via} Spearman's $\rho$ values")

    # add title to the figure
    clustermap.figure.suptitle(f"Biomarker Fold Change\n{title_str}", fontsize = "xx-large")

    plt.show()                                                      # show the image
    plt.close("all")
    return

def save_clustermap(sort_via: str, primary_only: bool = False, separate_biomarkers: bool = False, markers: str | None = None,
                    significant_biomarkers_only: bool = False, max_pval: float | None = None) -> None:
    """
    Saves a clustered heatmap from the dataset
    
    Parameters
    ----------
    sort_via: str
        metric to sort by
    primary_only: bool
        plot data from primary group only or primary and secondary groups
    separate_biomarkers: bool
        plot data from only one week of biomarkers or both weeks of biomarkers
    markers: str | None
        choose which week of biomarkers to plot if separate_biomarkers == True
        possible values: "4" or "20"
    significant_biomarkers_only: bool
        choose whether to only view biomarkers with a p-value below a desired value or not
    max_pval: float | None
        choose what is the maximum p-value to display if significant_biomarkers_only == True
        0 <= max_pval <= 1
    """

    save_dir = f"heatmap_change{markers}wk_{''.join(sort_via.split())}" # ex: heatmap_change20wk_TAT4wk

    # read in our datasets
    val_df = pd.read_csv(fold_change_path)                              # contains fold change data per patient and biomarker
    pval_df = pd.read_csv("modified_data/rho_pvals.csv")                # contains p-values of biomarkers per metric

    if(primary_only):
        val_df = val_df[val_df["in_primary"] == 1]                      # only use patients in the primary group by filtering out non-primary patients
        title_str = "Primary Group Only"
        save_dir += "_primary"
    else:
        title_str = "Primary and Secondary Groups"
        save_dir += "_both"


    # sort and extract relevant p-value data
    pval_df = pval_df[pval_df["exposure"] == sort_via]                  # find p-vals corresponding to the metric we want to sort by

    if(significant_biomarkers_only):                                    # only want to look at biomarkers that are significant
        pval_df = pval_df[pval_df["pval"] < max_pval]                   # keep rows with p-values below the level we choose
        biomarkers_to_keep = list(pval_df["biomarker"])                 # get list of the biomarkers that match this criteria
        patient_info_df = val_df.iloc[:, 0:6]
        val_df = val_df.loc[:, biomarkers_to_keep]                      # only keep the data from the biomarkers that matched criteria
        val_df = pd.concat([patient_info_df, val_df], axis = 1)

        save_dir = "figures/heatmap_significant/" + save_dir + f"_{max_pval}pval"

    else:
        save_dir = "figures/heatmap_all/" + save_dir

    pval_df = pval_df.sort_values(by = ["pval"])                        # sort p-vals from lowest to highest
    p_vals = pval_df["pval"]                                            # extract pd.Series of sorted p-values
    rhos = pval_df["rho"]                                               # extract pd.Series of rho values

    # keep biomarkers as determined by the flag if we want to separate them: 4 week or 20 week
    if(separate_biomarkers == True):
        patient_info_df = val_df.iloc[:, 0:6]                           # splitting up the dataframe into the patient info and the actual data
        val_df = val_df.iloc[:, 6:]
        markers = "." + markers                                         # minor string modification         #type: ignore
        cols_to_keep = val_df.columns.str.endswith(markers)             # creates boolean array corresponding to which columns are from the correct week    #type: ignore
        val_df = val_df.loc[:, cols_to_keep]                            # grabs dataframe with the correct columns only
        val_df = pd.concat([patient_info_df, val_df], axis = 1)

    # sort and extract relevant biomarker and fold change data
    num_patients, num_cols = val_df.shape                               # find number of patients (rows) and number of columns
    num_biomarkers = num_cols - 6                                       # num_biomarkers = num_cols - number of non-biomarker columns
    ids = val_df["Patient ID"]                                          # extract pd.Series of patient IDs (Patient ID column)
    col_to_sort = list(val_df.columns).index(sort_via)                  # find column that we want to use to sort by
    indices = [0, col_to_sort] + list(range(6, num_cols))               # get indices of columns to take the patient IDs, metric to sort by, and all data points
    val_df = val_df.iloc[:, indices]                                    # grab those columns in our data, toss out the rest
    val_df = val_df.fillna(0)                                           # replace NaNs (Not a Number) with 0 since sns.clustermap cannot have NaNs

    val_df = val_df.sort_values(by = [sort_via])                        # sort by sort_via metric: lowest -> highest
    exposure = val_df[sort_via]                                         # extract pd.Series of exposure metrics

    df_plot = val_df.transpose()                                        # transpose our dataframe: a column is now a patient
    df_plot.columns = df_plot.iloc[0]                                   # set our column labels to be our patient IDs
    df_plot = df_plot.iloc[2:]                                          # get rid of the patient IDs row and exposure metric rows

    # create objects for coloring patient exposure, p-values, and rho values
    exposure_cmap = sns.color_palette("flare", as_cmap = True)                                                      # create exposure colormap (cmap) object
    exposure_norm = matplotlib.colors.Normalize(vmin = exposure.min(), vmax = exposure.max())                       # create exposure normalizer 
    exposure_mapper = matplotlib.cm.ScalarMappable(norm = exposure_norm, cmap = exposure_cmap) #type: ignore        # create mapper which takes in a patient's exposure and outputs a color
    exposure_colors = [exposure_mapper.to_rgba(pat) for pat in exposure]                                            # create list of colors, from lowest exposure to highest

    # create color palette series: df_plot.columns contains sorted patient IDs, used as dict keys to index into the palette
    patient_color_palette = pd.Series(data = exposure_colors, index = df_plot.columns, name = f"{sort_via} Exposure")

    # create colors for p-values
    pval_cmap = sns.color_palette("viridis", as_cmap = True)
    if(significant_biomarkers_only):                                                                                # set vmax to maximum value allowed for p-values
        pval_norm = matplotlib.colors.Normalize(vmin = 0, vmax = max_pval)                                          # when looking at significant biomarkers only
    else:
        pval_norm = matplotlib.colors.Normalize(vmin = 0, vmax = 1)
    pval_mapper = matplotlib.cm.ScalarMappable(norm = pval_norm, cmap = pval_cmap)  #type: ignore
    pval_colors = [pval_mapper.to_rgba(p_val) for p_val in p_vals]

    # repeat above except for rho values
    rho_cmap = sns.color_palette("BrBG", as_cmap = True)
    rho_norm = matplotlib.colors.Normalize(vmin = -0.75, vmax = 0.75)
    rho_mapper = matplotlib.cm.ScalarMappable(norm = rho_norm, cmap = rho_cmap)     #type: ignore
    rho_colors = [rho_mapper.to_rgba(rho) for rho in rhos]

    # create a pd.DataFrame that contains the corresponding p-values and rho values for each biomarker
    row_color_palette = pd.DataFrame({"p-values": pval_colors,
                                    r"Spearman's $\rho$": rho_colors},
                                    index = pval_df["biomarker"])

    # plots our clustered heatmap, using our dataframe from above and clustering via rows ONLY
    # patients have increasing exposure from left to right: lowest -> highest exposure according to set metric
    # clustering: hierarchal clustering w/ standardized euclidean distance used
    # NOTE: scipy must be installed for sns.clustermap to work
    clustermap = sns.clustermap(df_plot, metric = "seuclidean", row_cluster = True, col_cluster = False, linewidth = 0.5,
                                row_colors = row_color_palette, col_colors = patient_color_palette, cmap = "plasma",
                                xticklabels = True, yticklabels = True, figsize=(20, 12), #vmin = 0, vmax = 1, center = 0.5, 
                                cbar_kws = {"label": "Scaled Fold Change", "orientation": "horizontal"},
                                cbar_pos = (0.6, 0.9, 0.3, 0.025))

    if(not significant_biomarkers_only):        # reduce font size of biomarkers so they don't overlap
        clustermap.ax_heatmap.set_yticklabels(clustermap.ax_heatmap.get_ymajorticklabels(), fontsize = 6)

    # plots colorbars for p-values and biomarkers onto the figure
    exposure_ax = clustermap.figure.add_axes([0.2, 0.9, 0.3, 0.025])                                # add a new axes onto the figure for the colorbar
    exposure_cbar = plt.colorbar(mappable = exposure_mapper, cax = exposure_ax,                     # create a new color bar using the mapper
                                 orientation = "horizontal", label = f"{sort_via} Exposure")        # corresponding to that colorbar

    pval_ax = clustermap.figure.add_axes([0.05, 0.85, 0.02, 0.1])
    pval_cbar = plt.colorbar(mappable = pval_mapper, cax = pval_ax,
                            orientation = "vertical", label = f"{sort_via} p-values")

    rho_ax = clustermap.figure.add_axes([0.1, 0.85, 0.02, 0.1])
    rho_cbar = plt.colorbar(mappable = rho_mapper, cax = rho_ax, orientation = "vertical",
                            label = rf"{sort_via} Spearman's $\rho$ values")

    # add title to the figure
    clustermap.figure.suptitle(f"Biomarker Fold Change\n{title_str}", fontsize = "xx-large")

    save_dir += ".png"
    clustermap.figure.savefig(save_dir, dpi=400)
    plt.close("all")
    return

def save_all_clustermaps(pval = 0.05):
    metrics = ["TAT 4 wk", "TAT 20 wk", "Cmax", "Cmax4"]            # possible metrics to use
    markers = ["4", "20"]                                           # week of metrics taken

    for metric in metrics:
        for marker in markers:
            save_clustermap(sort_via = metric, primary_only = False, separate_biomarkers = True, markers = marker,
                            significant_biomarkers_only = False, max_pval = None)
            save_clustermap(sort_via = metric, primary_only = True, separate_biomarkers = True, markers = marker,
                            significant_biomarkers_only = False, max_pval = None)
            save_clustermap(sort_via = metric, primary_only = False, separate_biomarkers = True, markers = marker,
                            significant_biomarkers_only = True, max_pval = pval)
            save_clustermap(sort_via = metric, primary_only = True, separate_biomarkers = True, markers = marker,
                            significant_biomarkers_only = True, max_pval = pval)
    
    return

save_all_clustermaps(pval=0.25)