#!/usr/bin/env python3
"""
Step 01: Environment Setup - Install dependencies and create directory structure
"""

import os
import sys
import subprocess
from pathlib import Path

# ✅ List of required Python packages
REQUIRED_PACKAGES = [
    "pandas>=1.5.0",
    "numpy<2",
    "scikit-learn>=1.3.0",
    "imbalanced-learn>=0.11.0",
    "torch>=2.0.0",
    "transformers>=4.30.0",
    "matplotlib>=3.7.0",
    "seaborn>=0.12.0",
    "requests>=2.31.0",
    "tqdm>=4.65.0",
    "joblib>=1.3.0",
    "rdkit-pypi>=2022.9.5"
]

# ✅ Folder structure for Git-tracked multimodal pipeline
PROJECT_DIRS = [
    "data/raw",                 # For raw clinical trial files (.csv)
    "data/processed",           # Cleaned/engineered outputs
    "models",                   # Trained model artifacts (.pkl, .pt)
    "results",                  # Evaluation reports, metrics
    "exports",                  # Exported CSVs or predictions
    "visualizations",           # Confusion matrices, SHAP plots etc
    "documentation",           # Project docs, mapping sheets
    "examples",                 # Sample notebooks or outputs
    "deployment_scripts",      # Prediction APIs or offline inference
    "tests",                   # Unit tests if added later
    "notebooks",               # Interactive exploration or demo notebooks
    "src/data_processing",     # Source code: data ingestion, cleaning
    "src/fingerprinting",      # Source code: SMILES, RDKit
    "src/embeddings",          # Source code: BERT embeddings
    "src/multimodal",          # Source code: feature integrator + models
    "src/utils"                # Shared utilities
]

def install_dependencies():
    print("📦 Checking/installing dependencies...")
    for pkg in REQUIRED_PACKAGES:
        pkg_name = pkg.split(">=")[0].replace("-", "_")  # handle e.g. "scikit-learn"
        try:
            __import__(pkg_name)
        except ImportError:
            print(f"⏳ Installing {pkg}...")
            subprocess.check_call([sys.executable, "-m", "pip", "install", pkg])
    print("✅ All dependencies verified or installed.")

def create_directories():
    print("📂 Creating project folders...")
    for d in PROJECT_DIRS:
        Path(d).mkdir(parents=True, exist_ok=True)
    print("✅ Folder structure created.")

if __name__ == "__main__":
    install_dependencies()
    create_directories()
    print("🚀 Environment setup completed.\nReady to begin data collection.")
