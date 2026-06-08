# SOD1 Dataset

This folder contains the current SOD1 ASO inhibitory-efficacy preparation workflow.

## Inputs

- `ASOptimizer_main.csv` is the raw SOD1 source file.
- `training/SOD1_main.csv` is a copy of the full starting training dataset.
- `training/SOD1_cleaned_v1.csv`, `SOD1_cleaned_v2.csv`, and `SOD1_cleaned_v3.csv` are successive cleaned versions used by the notebooks.

## Generated Outputs

- `training/03_chemistry strings into position-level annotations/` contains sequence-length validation, position-level sugar annotations, linkage annotations, and chemistry summary features.
- `training/04_inhibition_groups_and_leakage_safe_splits/` contains inhibition groups and leakage-safe split labels.
- `training/05_transcript_coordinate_mapping/` contains SOD1 transcript matching results and the fetched RefSeq FASTA.
- `training/06_biophysical_features/` contains sequence, gap, k-mer, transcript-coordinate, melting-temperature, and RDKit features.
- `training/07_model_ready_dataset/` contains the final model-ready dataset, its data dictionary, and validation summary.

## Code

The original notebooks are kept in `python/` for provenance. The cleaned, documented pipeline code is in `../../src/aso_prepare/`.
