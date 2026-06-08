"""Project-relative filesystem paths used by the ASO preparation code."""

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATASET_DIR = PROJECT_ROOT / "dataset"

SOD1_DIR = DATASET_DIR / "SOD1"
SOD1_TRAINING_DIR = SOD1_DIR / "training"
SOD1_VALIDATION_DIR = SOD1_DIR / "validation"

SOD1_RAW = SOD1_DIR / "ASOptimizer_main.csv"
SOD1_CLEANED_V1 = SOD1_TRAINING_DIR / "SOD1_cleaned_v1.csv"
SOD1_CLEANED_V2 = SOD1_TRAINING_DIR / "SOD1_cleaned_v2.csv"
SOD1_CLEANED_V3 = SOD1_TRAINING_DIR / "SOD1_cleaned_v3.csv"

SOD1_CHEMISTRY_DIR = SOD1_TRAINING_DIR / "03_chemistry strings into position-level annotations"
SOD1_SPLIT_DIR = SOD1_TRAINING_DIR / "04_inhibition_groups_and_leakage_safe_splits"
SOD1_TRANSCRIPT_DIR = SOD1_TRAINING_DIR / "05_transcript_coordinate_mapping"
SOD1_FEATURE_DIR = SOD1_TRAINING_DIR / "06_biophysical_features"
SOD1_MODEL_READY_DIR = SOD1_TRAINING_DIR / "07_model_ready_dataset"
