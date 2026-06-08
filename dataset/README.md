# Dataset Folder

This folder stores target-specific ASO datasets.

- `SOD1/` contains the current prepared SOD1 dataset, intermediate processing outputs, validation data, and original notebooks.
- `SMN2/` is reserved for future SMN2 preparation work. Use the same staged layout as SOD1 when SMN2 source data becomes available.

Generated CSVs are kept because they document the current data provenance. New reusable code should live under `src/aso_prepare/` rather than inside notebook checkpoint folders.
