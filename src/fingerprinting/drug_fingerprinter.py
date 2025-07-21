#!/usr/bin/env python3
"""
Step 03: Drug Fingerprint Generation
- Fetch SMILES from PubChem
- Generate Morgan fingerprints using RDKit
- Output to data/processed/unique_drugs_with_fingerprint.csv
"""

import pandas as pd
import requests
from rdkit import Chem
from rdkit.Chem import AllChem
from pathlib import Path
import time

def get_smiles_from_pubchem(drug_name):
    """Fetch canonical SMILES for a given drug name from PubChem."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/property/SMILES/JSON"
    try:
        response = requests.get(url, timeout=10)
        if response.ok:
            data = response.json()
            return data['PropertyTable']['Properties'][0]['SMILES']
    except Exception:
        pass
    return None

def get_fingerprint_from_smiles(smiles):
    """Generate RDKit Morgan fingerprint as a string for a SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        return fp.ToBitString()
    return None

def main(
    drugs_file="data/processed/unique_drugs.csv",
    output_file="data/processed/unique_drugs_with_fingerprint.csv"
):
    drugs_df = pd.read_csv(drugs_file)
    smiles_list = []
    fingerprint_list = []

    for idx, row in drugs_df.iterrows():
        name = row["drug_name"]
        smiles = get_smiles_from_pubchem(name)
        if smiles:
            fingerprint = get_fingerprint_from_smiles(smiles)
        else:
            fingerprint = None
        smiles_list.append(smiles)
        fingerprint_list.append(fingerprint)
        time.sleep(0.2)  # Polite delay to avoid rate-limiting

    drugs_df["smiles"] = smiles_list
    drugs_df["fingerprint"] = fingerprint_list
    drugs_df = drugs_df[drugs_df["fingerprint"].notnull()]  # Keep only those with valid FP
    drugs_df.to_csv(output_file, index=False)
    print(f"✅ Saved {len(drugs_df)} drugs with SMILES and fingerprints to {output_file}")

if __name__ == "__main__":
    main()
