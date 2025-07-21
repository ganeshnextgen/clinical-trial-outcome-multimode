#!/usr/bin/env python3
"""
Step 01: Environment Setup - Installs dependencies and prepares directory structure.
"""

import os
import sys
import subprocess
from pathlib import Path

# List of required packages
REQUIRED_PACKAGES = [
    "pandas>=1.5.0",
    "numpy>=1.24.0",
    "scikit-learn>=1.3.0",
    "imbalanced-learn>=0.11.0",
    "torch>=2.0.0",
    "transformers>=4.30.0",
    "matplotlib>=3.7.0",
    "seaborn>=0.12.0",
    "requests>=2.31.0",
    "tqdm>=4.65.0",
    "joblib>=1.3.0"
]

def install_and_confirm():
    print("📦 Installing required packages...")
    for pkg in REQUIRED_PACKAGES:
        pkg_name = pkg.split(">=")[0]
        try:
            __import__(pkg_name.replace("-", "_"))
        except ImportError:
            subprocess.check_call([sys.executable, "-m", "pip", "install", pkg, "-q"])
    print("✅ All dependencies installed.")

def create_project_dirs():
    dirs = [
        "data/processed",
        "models",
        "results",
        "exports",
        "visualizations",
        "documentation",
        "examples",
        "deployment_scripts",
        "tests"
    ]
    for d in dirs:
        Path(d).mkdir(parents=True, exist_ok=True)
    print("✅ Directory structure created.")

if __name__ == "__main__":
    install_and_confirm()
    create_project_dirs()
    print("🚀 Environment setup completed.")
