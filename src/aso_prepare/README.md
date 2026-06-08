# `aso_prepare` Code Guide

This package turns raw or cleaned ASO assay tables into model-ready feature tables.

## Files

- `paths.py` centralizes project-relative paths. This replaces hard-coded local Windows paths in the original notebooks.
- `sod1_pipeline.py` contains the SOD1 preparation workflow as documented functions. Each function corresponds to one logical notebook stage and explains its inputs, outputs, and assumptions.
- `__init__.py` exposes the package version.

## Main Concepts

- `Sequence` is normalized to uppercase DNA alphabet characters.
- `Chemical_Pattern` is interpreted position by position: `M` means MOE, `C` means cEt, and `d` means DNA.
- Linkage locations are treated as 1-based internucleotide positions.
- Sequence leakage is controlled by grouping train/test splits by `Sequence`.
- Optional dependencies are handled gracefully:
  - Biopython is used for transcript retrieval and DNA melting-temperature features.
  - RDKit is used for SMILES-derived molecular descriptors.

## Example

From the repository root:

```bash
python -m src.aso_prepare.sod1_pipeline
```

The command reruns the SOD1 feature-building pipeline using project-relative paths and writes outputs under `dataset/SOD1/training/`.
