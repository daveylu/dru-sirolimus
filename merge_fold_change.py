import numpy as np
import pandas as pd

def relabel_4wk(col_label):
    if("foldchange" in col_label):
        return col_label + ".4"
    else:
        return col_label

def relabel_20wk(col_label):
    if("foldchange" in col_label):
        return col_label + ".20"
    else:
        return col_label

cd_path = "./modified_data/cd4cd8_data.csv"
cd_df = pd.read_csv(cd_path)

nuc_path = "./modified_data/rna_dna_data.csv"
nuc_df = pd.read_csv(nuc_path)

solubles_path = "./modified_data/soluble_biomarker_data.csv"
solubles_df = pd.read_csv(solubles_path)

unstim_path = "./modified_data/unstim_flow_data.csv"
unstim_df = pd.read_csv(unstim_path)

stim_path = "./modified_data/stim_flow_data.csv"
stim_df = pd.read_csv(stim_path)

df = pd.merge(nuc_df, solubles_df, on = ["patid", "rxcalwk"])
df = pd.merge(df, cd_df, on = ["patid", "rxcalwk"])
df = pd.merge(df, unstim_df, on = ["patid", "rxcalwk"])
df = pd.merge(df, stim_df, on = ["patid", "rxcalwk"])

df = df[df["rxcalwk"].isin([4, 20])]
df_front = df.iloc[:, 0:2]
df_foldchange = df.iloc[:, df.columns.str.contains("foldchange")]
df = pd.concat([df_front, df_foldchange], axis = 1)

markers = df_foldchange.columns # type: ignore
df_4wk = df[df["rxcalwk"] == 4]
df_4wk = df_4wk.drop("rxcalwk", axis = 1)
df_20wk = df[df["rxcalwk"] == 20]
df_20wk = df_20wk.drop("rxcalwk", axis = 1)

df_4wk = df_4wk.rename(mapper = relabel_4wk, axis = 1)
df_20wk = df_20wk.rename(mapper = relabel_20wk, axis = 1)

df = pd.merge(df_4wk, df_20wk, on = "patid")

df.to_csv("./modified_data/fold_change.csv")
