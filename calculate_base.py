import numpy as np
import pandas as pd

def calculate_base(df, save_path):

    patids = df["patid"].unique()
    biomarkers = list(df.columns)[2:]

    for biomarker in biomarkers:
        base_vals = []
        for patid in patids:
            pat_df = df[df["patid"] == patid]
            num_to_append = pat_df.shape[0]
            pat_df = pat_df[pat_df["rxcalwk"].isin([-12, 0])]
            biomarker_vals = pat_df[biomarker].to_numpy()
            base_val = np.round(np.mean(biomarker_vals), 6)
            base_vals += [base_val] * num_to_append
        df[f"base_{biomarker}"] = base_vals

    df.to_csv(save_path, index=False)

cd_path = "./modified_data/cd4cd8_data.csv"
cd_df = pd.read_csv(cd_path)
calculate_base(cd_df, cd_path)

nuc_path = "./modified_data/rna_dna_data.csv"
nuc_df = pd.read_csv(nuc_path)
calculate_base(nuc_df, nuc_path)

solubles_path = "./modified_data/soluble_biomarker_data.csv"
solubles_df = pd.read_csv(solubles_path)
calculate_base(solubles_df, solubles_path)

unstim_path = "./modified_data/unstim_flow_data.csv"
unstim_df = pd.read_csv(unstim_path)
calculate_base(unstim_df, unstim_path)

stim_path = "./modified_data/stim_flow_data.csv"
stim_df = pd.read_csv(stim_path)
calculate_base(stim_df, stim_path)