#!/usr/bin/env python3
"""
Step 03: Data Preprocessing & Feature Engineering
"""

import pandas as pd
import numpy as np
import re
from sklearn.preprocessing import LabelEncoder, StandardScaler
from pathlib import Path

class ClinicalTrialPreprocessor:
    def __init__(self, data_dir="data", processed_dir="data/processed"):
        self.data_dir = Path(data_dir)
        self.processed_dir = Path(processed_dir)
        self.processed_dir.mkdir(parents=True, exist_ok=True)
        self.label_encoders = {}

    def load_raw_files(self):
        files = sorted(self.data_dir.glob("raw_clinical_trials_*.csv"))
        if not files:
            raise FileNotFoundError("No raw files found.")
        
        all_dfs = []
        for file in files:
            print(f"📂 Loading: {file}")
            df = pd.read_csv(file)
            all_dfs.append(df)
        
        combined_df = pd.concat(all_dfs, ignore_index=True)
        print(f"✅ Loaded and combined {len(files)} files.")
        return combined_df


    def preprocess(self, df):
        print("🧹 Cleaning and encoding data...")
    
        # Case-insensitive matching for outcome
        succ = ['completed', 'active, not recruiting', 'recruiting', 'enrolling by invitation']
        fail = ['terminated', 'withdrawn', 'suspended', 'no longer available', 'temporarily not available']
    
        def outcome(status, reason):
            if pd.isna(status):
                return np.nan
            s = str(status).strip().lower()
            if s in succ:
                return 1
            if s in fail:
                return 0
            if pd.notna(reason) and str(reason).strip() != "":
                return 0
            return np.nan
    
        df['outcome'] = df.apply(lambda r: outcome(r.get('overall_status'), r.get('why_stopped')), axis=1)
        df = df.dropna(subset=['outcome']).copy()
        
        print(f"✅ Outcome labeling complete. Remaining trials: {len(df)}")

        # Text cleaning and join
        text_cols = ['brief_title', 'official_title', 'brief_summary', 'detailed_description']
        def clean(x):
            if pd.isnull(x): return ""
            return re.sub(r"[^\w\s]", " ", str(x).lower())
        df['combined_text'] = df[text_cols].fillna("").agg(' '.join, axis=1).apply(clean)
        df['enrollment_count'] = pd.to_numeric(df['enrollment_count'], errors='coerce').fillna(0)

        # Encode categoricals
        cat_cols = ['study_type', 'phases', 'allocation', 'intervention_model', 'masking', 'primary_purpose']
        for col in cat_cols:
            le = LabelEncoder()
            df[f"{col}_encoded"] = le.fit_transform(df[col].fillna("Unknown"))
            self.label_encoders[col] = le

        # Save processed
        out_file = self.processed_dir / "clinical_trials_processed.csv"
        df.to_csv(out_file, index=False)
        print(f"✅ Saved processed data: {out_file}")
        return df

if __name__ == "__main__":
    processor = ClinicalTrialPreprocessor()
    raw = processor.load_raw_files()
    _ = processor.preprocess(raw)
    print("🚀 Preprocessing completed.")
