import pandas as pd
import numpy as np

def fold_change(df, save_path):
    biomarkers = []
    col_names = list(df.columns)[2:]
    for name in col_names:
        if("base" not in name):
            biomarkers.append(name)

    d = dict()

    for biomarker in biomarkers:
        biomarker = biomarker
        val = df[biomarker].to_numpy()
        base = df["base_" + biomarker].to_numpy()
        fold_change = np.nan_to_num(val / base, nan=1, posinf=1, neginf=1)
        d[biomarker + "_foldchange"] = fold_change

    temp = pd.DataFrame(d)

    df = pd.concat((df, temp), axis = 1)
    df.to_csv(save_path, index=False)

cd_path = "./modified_data/cd4cd8_data.csv"
cd_df = pd.read_csv(cd_path)
fold_change(cd_df, cd_path)

nuc_path = "./modified_data/rna_dna_data.csv"
nuc_df = pd.read_csv(nuc_path)
fold_change(nuc_df, nuc_path)

solubles_path = "./modified_data/soluble_biomarker_data.csv"
solubles_df = pd.read_csv(solubles_path)
fold_change(solubles_df, solubles_path)

unstim_path = "./modified_data/unstim_flow_data.csv"
unstim_df = pd.read_csv(unstim_path)
fold_change(unstim_df, unstim_path)

stim_path = "./modified_data/stim_flow_data.csv"
stim_df = pd.read_csv(stim_path)
fold_change(stim_df, stim_path)