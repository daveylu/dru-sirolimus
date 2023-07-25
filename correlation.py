import numpy as np
import pandas as pd
import scipy.stats

df = pd.read_csv("./modified_data/fold_change.csv")

exposures = ["TAT 4 wk", "TAT 20 wk", "Cmax", "Cmax4"]
biomarkers = list(df.columns)[6:]

d = dict()
save = pd.DataFrame()

for metric in exposures:
    exposure = df[metric]
    metric_list = [metric] * len(biomarkers)
    for biomarker in biomarkers:
        marker = df[biomarker]
        result = scipy.stats.spearmanr(a = exposure, b = marker, nan_policy = "omit") # type: ignore
        d[biomarker] = result

    temp = pd.DataFrame(d).T
    temp.columns = ["rho", "pval"]
    temp = temp.reset_index(names = "biomarker")
    temp.insert(0, "exposure", metric_list)

    save = pd.concat([save, temp], axis = 0)


save.to_excel("./modified_data/temp.xlsx")
