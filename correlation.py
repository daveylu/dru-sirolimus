import numpy as np
import pandas as pd
import scipy.stats

df = pd.read_csv("./modified_data/fold_change.csv")

exposures = ["TAT 4 wk", "TAT 20 wk", "Cmax", "Cmax4"]
biomarkers = list(df.columns)[6:]

d = dict()
save = pd.DataFrame()

for metric in exposures:
    exposure_both = df[metric]
    exposure_primary = df[metric][df["in_primary"] == 1]
    metric_list = [metric] * len(biomarkers)
    for biomarker in biomarkers:
        # primary and secondary
        marker_both = df[biomarker]
        result_both = scipy.stats.spearmanr(a = exposure_both, b = marker_both, nan_policy = "omit") # type: ignore

        # primary only
        marker_primary = df[biomarker][df["in_primary"] == 1]
        result_primary = scipy.stats.spearmanr(a = exposure_primary, b = marker_primary, nan_policy = "omit") # type: ignore
        d[biomarker] = result_both + result_primary

    temp = pd.DataFrame(d).T
    temp.columns = ["rho_both", "pval_both", "rho_primary", "pval_primary"]
    temp = temp.reset_index(names = "biomarker")
    temp.insert(0, "exposure", metric_list)

    save = pd.concat([save, temp], axis = 0)

save.to_csv("./modified_data/rho_pvals.csv", index=False)
