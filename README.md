# ASO Inhibitory Efficacy Framework

This repository prepares machine-learning datasets for Antisense Oligonucleotides (ASOs) targeting SOD1 and SMN2. The current implemented pipeline is for SOD1. SMN2 is present as a target folder placeholder, but does not yet contain source data or preprocessing code.

## Current Structure

- `dataset/SOD1/ASOptimizer_main.csv` - raw multi-gene ASO assay export (31,752 rows across 21 target genes) used as the starting dataset; the SOD1 subset is selected from it during cleaning.
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

## Known Limitations

Please read these before training models on the current datasets:

- **SMILES do not encode the actual ASO chemistry.** The committed `Smiles` strings represent plain phosphodiester DNA oligos: they contain no phosphorothioate (`P(=S)`) groups and no MOE/cEt sugar modifications, even for fully phosphorothioated, MOE/cEt-modified ASOs. Consequently, all `rdkit_*` descriptor columns are **placeholders describing unmodified DNA**, not the real molecules, and should not be used for modeling until chemistry-aware structures are generated.
- **5-methylcytosine is not encoded.** The `Modification` column indicates 5mC for all rows, but neither `Sequence` nor `Chemical_Pattern` carries this information; it is silently dropped from features.
- **Melting-temperature features are plain-DNA approximations.** `basic_DNA_Tm_NN_*` uses Biopython's DNA nearest-neighbor table on heavily modified oligos; treat these values as rough proxies only.
- **Train/test splitting prevents only exact-sequence leakage.** `split_2way` groups by the exact `Sequence` string. ASOs targeting overlapping transcript windows (shifted by 1-2 nt) can still appear on both sides of the split, so held-out metrics may be optimistically biased. A homology/similarity-aware split is recommended for future modeling.
- **The validation set overlaps the training data.** Several entries in `dataset/SOD1/validation/SOD1_validation_dataset_smiles.csv` share ISIS identifiers and/or exact sequences with rows in the training table. Do not treat it as an independent external validation set in its current form.
- **The packaged pipeline does not yet regenerate the committed CSVs end-to-end** (the SOD1 gene filter and the NCBI transcript-mapping stage are not yet ported from the exploratory notebooks). The notebooks under `dataset/SOD1/python/` are kept for provenance and contain hard-coded local paths; use `src/aso_prepare/` for reruns.

## Dataset Provenance

See `dataset/PROVENANCE.md` for the (in-progress) record of the source, citation, and filtering rules for the raw export.
