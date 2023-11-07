import pandas as pd
from sklearn.preprocessing import minmax_scale

df = pd.read_csv("modified_data/fold_change_no_outliers.csv")
front = df.iloc[:, :6]
df = df.iloc[:, 6:]
biomarkers = list(df.columns)

a = df.to_numpy()
a = minmax_scale(a, axis = 0)   # scales data to a [0, 1] range

df = pd.DataFrame(a, columns = biomarkers)

df = pd.concat([front, df], axis = 1)

df.to_csv("modified_data/scaled_fold_change_no_outliers.csv", index = False)