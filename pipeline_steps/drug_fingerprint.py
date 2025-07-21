# drug_fingerprint.py

import pandas as pd
import requests
from rdkit import Chem
from rdkit.Chem import AllChem
from pathlib import Path

def get_smiles_from_pubchem(drug_name):
    """Fetch canonical SMILES for a given drug name from PubChem."""
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/property/SMILES/JSON'
    response = requests.get(url)
    if response.ok:
        try:
            data = response.json()
            return data['PropertyTable']['Properties'][0]['SMILES']
        except (KeyError, IndexError):
            return None
    return None

def get_fingerprint_from_smiles(smiles):
    """Generate RDKit Morgan fingerprint as a string for a SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        return fp.ToBitString()
    return None

def main(drugs_file="data/processed/unique_drugs.csv", output_file="data/processed/unique_drugs_with_fingerprint.csv"):
    drugs_df = pd.read_csv(drugs_file)
    drugs_df['smiles'] = drugs_df['drug_name'].apply(get_smiles_from_pubchem)
    drugs_df['fingerprint'] = drugs_df['smiles'].apply(
        lambda s: get_fingerprint_from_smiles(s) if pd.notnull(s) else None
    )
    # Optionally keep only those with a fingerprint
    drugs_df = drugs_df[drugs_df['fingerprint'].notnull()]
    drugs_df.to_csv(output_file, index=False)
    print(f"✅ Saved {len(drugs_df)} drugs with SMILES and fingerprint.")

if __name__ == "__main__":
    main()
