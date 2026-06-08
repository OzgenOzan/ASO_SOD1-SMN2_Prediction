# ASO Inhibitory Efficacy Framework

This repository prepares machine-learning datasets for Antisense Oligonucleotides (ASOs) targeting SOD1 and SMN2. The current implemented pipeline is for SOD1. SMN2 is present as a target folder placeholder, but does not yet contain source data or preprocessing code.

## Current Structure

- `dataset/SOD1/ASOptimizer_main.csv` - raw SOD1 export used as the starting dataset.
- `dataset/SOD1/training/` - cleaned and feature-engineered SOD1 training datasets.
- `dataset/SOD1/validation/` - external or held-out validation data.
- `dataset/SOD1/python/` - original exploratory notebooks used to create the current CSV outputs.
- `src/aso_prepare/` - cleaned, documented Python code for rerunning the SOD1 preparation workflow.
- `dataset/SMN2/` - placeholder for the future SMN2 dataset preparation workflow.

## SOD1 Pipeline

The cleaned code follows the existing notebook-derived stages:

1. Select one assay record per ASO/location/SMILES combination, keeping the highest ASO dose and longest treatment duration.
2. Validate sequence length and chemistry-pattern length.
3. Convert chemistry strings into position-level and linkage-level annotations.
4. Create inhibition bins and sequence-leakage-safe train/test splits.
5. Add transcript-coordinate mapping outputs when available.
6. Add sequence, gap, k-mer, transcript-coordinate, melting-temperature, and RDKit descriptor features.
7. Export the model-ready SOD1 dataset and validation summary.

See `src/aso_prepare/README.md` for code-level explanations.

## Notes

The current final SOD1 model-ready file is:

`dataset/SOD1/training/07_model_ready_dataset/SOD1_model_ready_v1.csv`

It contains 1,026 rows and 162 columns. The target variable is `Inhibition(%)`.
