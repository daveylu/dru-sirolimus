import pandas as pd
import numpy as np

df = pd.read_csv("modified_data/fold_change_no_outliers.csv")
front = df.iloc[:, :6]
df = df.iloc[:, 6:]

df = np.log2(df)

col_names = list(df.columns)
new_names = []

for col_name in col_names:
    index = col_name.find("foldchange")
    new_name = col_name[:index] + "log" + col_name[index:]
    new_names.append(new_name)

df.columns = new_names
df = pd.concat([front, df], axis = 1)

df.to_csv("modified_data/log_fold_change_no_outliers.csv", index = False)