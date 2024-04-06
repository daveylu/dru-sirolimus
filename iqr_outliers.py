import numpy as np
import pandas as pd

df = pd.read_csv("modified_data/fold_change_no_tinyvals.csv")
head = df.iloc[:, 0:6]
df = df.iloc[:, 6:]
col_names = list(df.columns)

a = df.to_numpy()
min = np.percentile(a, q = 2.5, axis = 0)
max = np.percentile(a, q = 97.5, axis = 0)

low_outliers = a < min
high_outliers = a > max

rows, cols = a.shape

for i in range(rows):
    for j in range(cols):
        if(high_outliers[i, j] == True):
            not_high_outliers = a[:, j][np.invert(high_outliers[:, j])]
            a[i, j] = np.max(not_high_outliers)
        elif(low_outliers[i, j] == True):
            not_low_outliers = a[:, j][np.invert(low_outliers[:, j])]
            a[i, j] = np.min(not_low_outliers)

new_df = pd.DataFrame(a, columns = col_names)
df = pd.concat([head, new_df], axis = 1)
df.to_csv("modified_data/fold_change_no_outliers.csv", index = False)