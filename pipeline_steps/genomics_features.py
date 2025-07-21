#!/usr/bin/env python3
"""
genomics_features.py: Reduce and align genomics vectors for fusion
"""
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA

def process_genomics(input_csv="data/genomics_data.csv", output_file="data/processed/genomics_features.npy", n_components=50):
    df = pd.read_csv(input_csv)
    gene_features = df.filter(regex="^gene_").values
    pca = PCA(n_components=n_components, random_state=42)
    reduced = pca.fit_transform(gene_features)
    np.save(output_file, reduced)
    trial_ids = df[['trial_id']].copy()
    trial_ids.to_csv(output_file.replace(".npy", "_ids.csv"), index=False)
    print(f"Genomics features saved: {output_file} (PCA {n_components}D)")
    print(f"Explained variance ratio: {pca.explained_variance_ratio_.sum():.2f}")

if __name__ == "__main__":
    process_genomics()
