import pandas as pd
from sklearn.preprocessing import minmax_scale, robust_scale

df = pd.read_csv("modified_data/fold_change.csv")
front = df.iloc[:, :6]
df = df.iloc[:, 6:]
biomarkers = list(df.columns)

a = df.to_numpy()
a = robust_scale(a, axis = 0)

df = pd.DataFrame(a, columns = biomarkers)

df = pd.concat([front, df], axis = 1)

df.to_csv("modified_data/temp.csv", index = False)