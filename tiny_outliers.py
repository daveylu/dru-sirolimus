import numpy as np
import pandas as pd

def too_small(x):
    if(x < 0.05):
        return 1.0
    else:
        return x

df = pd.read_csv("modified_data/fold_change.csv")
head = df.iloc[:, 0:6]
df = df.iloc[:, 6:]
df = df.map(too_small)

df = pd.concat([head, df], axis = 1)
df.to_csv("modified_data/fold_change_no_tinyvals.csv", index=False)