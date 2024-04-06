import numpy as np
import pandas as pd
from scipy.stats import ttest_rel

data = pd.read_csv("./modified_data/data.csv")

# getting biomarker names...
col_names = np.array(data.columns)

is_base = []
for name in col_names:
    if("base" in name):
        is_base.append(True)
    else:
        is_base.append(False)
is_base = np.array(is_base)

base_names = col_names[is_base]
biomarkers = []
for base_name in base_names:
    biomarker = base_name[len("base_"):]
    biomarkers.append(biomarker)

names = []
t_stat = []
p_val = []
deg_free = []

for biomarker in biomarkers:
    # obtaining the corresponding labels for this biomarker
    base = "base_" + biomarker
    wk4 = biomarker + ".4"
    wk20 = biomarker + ".20"

    #obtaining the data from those labels
    base_data = data[base].to_numpy()
    wk4_data = data[wk4].to_numpy()
    wk20_data = data[wk20].to_numpy()

    # performing the t-tests
    wk4_result = ttest_rel(base_data, wk4_data)
    wk20_result = ttest_rel(base_data, wk20_data)

    # appending results to our lists
    names.append(wk4)
    t_stat.append(wk4_result.statistic)
    p_val.append(wk4_result.pvalue)
    deg_free.append(wk4_result.df)

    names.append(wk20)
    t_stat.append(wk20_result.statistic)
    p_val.append(wk20_result.pvalue)
    deg_free.append(wk20_result.df)

df = pd.DataFrame({
    "biomarker": names,
    "p-value": p_val,
    "t-statistic": t_stat,
    "degrees of freedom": deg_free
})

df.to_csv("./modified_data/t-tests.csv")
