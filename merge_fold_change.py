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

rna_dna = "./modified_data/rna_dna_data.xlsx"
solubles = "./modified_data/soluble_biomarker_data.xlsx"
unstim = "./modified_data/unstim_flow_data.xlsx"
cd = "./modified_data/cd4cd8_data.xlsx"
stim = "./modified_data/stim_flow_data.xlsx"

nuc_df = pd.read_excel(rna_dna, sheet_name = "Data")
sol_df = pd.read_excel(solubles, sheet_name = "Data")
cd_df = pd.read_excel(cd, sheet_name = "Data")
unstim_df = pd.read_excel(unstim, sheet_name = "Data")
stim_df = pd.read_excel(stim, sheet_name = "Data")

df = pd.merge(nuc_df, sol_df, on = ["patid", "rxcalwk"])
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

df.to_excel("./modified_data/temp.xlsx")
