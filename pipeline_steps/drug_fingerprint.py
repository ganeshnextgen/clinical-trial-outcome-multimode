#!/usr/bin/env python3
"""
Drug Fingerprint Generator for Clinical Trials (supports integration with pipeline)
"""

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys

class DrugFingerprintGenerator:
    def __init__(self, smiles_lookup_dict=None, fp_type="morgan", fp_size=2048, radius=2):
        """
        smiles_lookup_dict: dict mapping lowercase drug names to SMILES strings.
        fp_type: "morgan" (ECFP-like) or "maccs"
        fp_size: size for Morgan/ECFP fingerprints
        radius: radius for Morgan fingerprints
        """
        self.smiles_lookup = smiles_lookup_dict or {}
        self.fp_type = fp_type
        self.fp_size = fp_size
        self.radius = radius

    def fetch_smiles(self, drug_name):
        """
        Lookup or fetch SMILES. For deployment, use external API (e.g., PubChem) if not found.
        """
        return self.smiles_lookup.get(drug_name.lower(), "")

    def compute_fingerprint(self, smiles):
        """
        Returns a numpy binary vector for the fingerprint or None if invalid.
        """
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None
        if self.fp_type == "morgan":
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, self.radius, nBits=self.fp_size)
        elif self.fp_type == "maccs":
            fp = MACCSkeys.GenMACCSKeys(mol)
            return np.array(list(fp.ToBitString()[1:]), dtype=np.int8)  # skip first empty bit
        else:
            raise ValueError("Unsupported fingerprint type")
        arr = np.zeros((1,), dtype=np.int8)
        DataStructs.ConvertToNumpyArray(fp, arr)
        return arr

    def process_drug_list(self, drug_df):
        """
        Takes a DataFrame with 'drug_name' and optional 'smiles' columns.
        Returns drug_name, SMILES, and fingerprint features.
        """
        records = []
        for ix, row in drug_df.iterrows():
            drug = row['drug_name']
            smiles = row.get('smiles', '') or self.fetch_smiles(drug)
            if not smiles:
                continue
            fp = self.compute_fingerprint(smiles)
            if fp is not None:
                rec = {'drug_name': drug, 'smiles': smiles}
                for i, v in enumerate(fp):
                    rec[f"fp_{i}"] = v
                records.append(rec)
        return pd.DataFrame(records)

if __name__ == "__main__":
    # Example: Load drug names from preprocessing output
    drug_df = pd.read_csv("data/processed/unique_drugs.csv")
    # If you have {drug_name: smiles} mapping, supply it here
    # For demo, let's fill all SMILES as blank (user should update)
    drug_df['smiles'] = ""  # Fill with your own SMILES mapping if available

    gen = DrugFingerprintGenerator(fp_type="morgan", fp_size=1024)
    fp_df = gen.process_drug_list(drug_df)
    fp_df.to_csv("data/processed/drug_fingerprints.csv", index=False)
    print(f"✅ Generated and saved {len(fp_df)} drug fingerprints.")
