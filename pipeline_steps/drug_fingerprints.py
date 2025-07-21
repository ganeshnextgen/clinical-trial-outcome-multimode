#!/usr/bin/env python3
"""
drug_fingerprints.py: Encode SMILES into molecular vectors (e.g., ECFP with RDKit)
"""
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from pathlib import Path

def smiles_to_ecfp(smiles, radius=2, n_bits=2048):
    """Convert a SMILES string to ECFP fingerprint"""
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
        return np.array(fp)
    else:
        return np.zeros(n_bits)

def process_drugs(input_csv="data/drug_info.csv", output_file="data/processed/drug_fingerprints.npy"):
    df = pd.read_csv(input_csv)
    print(f"Loaded {len(df)} drugs with SMILES.")
    fps = np.stack([smiles_to_ecfp(s) for s in df['smiles']])
    np.save(output_file, fps)
    df_ids = df[['trial_id', 'drug_name']].copy()
    df_ids.to_csv(output_file.replace(".npy", "_ids.csv"), index=False)
    print(f"Drug fingerprints saved: {output_file}")

if __name__ == "__main__":
    process_drugs()
