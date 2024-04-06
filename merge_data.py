import numpy as np
import pandas as pd

cd_df = pd.read_csv("./modified_data/cd4cd8_data.csv")
nuc_df = pd.read_csv("./modified_data/rna_dna_data.csv")
solubles_df = pd.read_csv("./modified_data/soluble_biomarker_data.csv")
unstim_df = pd.read_csv("./modified_data/unstim_flow_data.csv")
stim_df = pd.read_csv("./modified_data/stim_flow_data.csv")

df = pd.merge(cd_df, nuc_df, how="inner", on=["patid", "rxcalwk"])
df = pd.merge(df, solubles_df, how="inner", on=["patid", "rxcalwk"])
df = pd.merge(df, unstim_df, how="inner", on=["patid", "rxcalwk"])
df = pd.merge(df, stim_df, how="inner", on=["patid", "rxcalwk"])

df_4wk = df[df["rxcalwk"] == 4]
df_20wk = df[df["rxcalwk"] == 20]

is_base = []
for col_name in list(df.columns):
    if("base" in col_name):
        is_base.append(True)
    else:
        is_base.append(False)

is_base = np.array(is_base)
is_not_base = ~is_base
is_base[0] = True           # this is just to keep the patid in the df
                            # so we can use it as a key to merge on

df_base = df_4wk.loc[:, is_base]
df_4wk = df_4wk.loc[:, is_not_base]
df_20wk = df_20wk.loc[:, is_not_base]

extra, biomarkers = np.array((df_4wk.columns))[:2], np.array((df_4wk.columns))[2:]
biomarkers_4wk = biomarkers + ".4"
biomarkers_20wk = biomarkers + ".20"

columns_4wk = np.concatenate((extra, biomarkers_4wk))
columns_20wk = np.concatenate((extra, biomarkers_20wk))

df_4wk.columns = columns_4wk
df_20wk.columns = columns_20wk

df_4wk = df_4wk.drop("rxcalwk", axis=1)
df_20wk = df_20wk.drop("rxcalwk", axis=1)

is_foldchange = []
for col_name in list(df_4wk.columns):
    if("foldchange" in col_name):
        is_foldchange.append(True)
    else:
        is_foldchange.append(False)

is_foldchange = np.array(is_foldchange)
is_not_foldchange = ~is_foldchange
is_foldchange[0] = True     # this is just to keep the patid in the df
                            # so we can use it as a key to merge on

df_4wk_change = df_4wk.loc[:, is_not_foldchange]
df_20wk_change = df_20wk.loc[:, is_not_foldchange]
df_4wk_foldchange = df_4wk.loc[:, is_foldchange]
df_20wk_foldchange = df_20wk.loc[:, is_foldchange]

df = pd.merge(df_base, df_4wk_change, how="inner", on="patid")
df = pd.merge(df, df_20wk_change, how="inner", on="patid")
df = pd.merge(df, df_4wk_foldchange, how="inner", on="patid")
df = pd.merge(df, df_20wk_foldchange, how="inner", on="patid")

df.to_csv("modified_data/data.csv", index=False)