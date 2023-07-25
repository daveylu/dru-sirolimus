import pandas as pd
import numpy as np

df = pd.read_excel("modified_data/stim_flow_data.xlsx", sheet_name = "Data")

biomarkers = list(df.columns)[2:38]

d = dict()

for biomarker in biomarkers:
    val = df[biomarker].to_numpy()
    base = df["base_" + biomarker].to_numpy()
    fold_change = np.nan_to_num(val / base, posinf = 0)
    d[biomarker + "_foldchange"] = fold_change

temp = pd.DataFrame(d)

df = pd.concat((df, temp), axis = 1)