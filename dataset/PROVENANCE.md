# Dataset Provenance

> **Status: TEMPLATE — TODO fields must be filled in by the data owner before this
> dataset is cited or relied upon.** This template was added during a repository
> audit; the audit could not recover the information below from the repository
> contents.

## Source

- Origin of `SOD1/ASOptimizer_main.csv` (database / publication / internal export): **TODO**
- Citation (paper, database entry, or URL): **TODO**
- Download / export date: **TODO**
- Version or accession of the source (if applicable): **TODO**
- License / terms of use of the underlying data: **TODO**

## Scope of the raw export

The raw file contains 31,752 assay rows spanning 21 target genes
(IRF4, HSD17B13, MALAT1, SOD-1, APOL1, DGAT2, SNHG14, SNCA, HBV, Yap1, IRF5,
PKK, Tau, K-RAS, SNCA_LNA, MYH7, ANGPTL2_LNA, HTRA1, HIF1A, plus controls).
Raw `Inhibition(%)` values range from -786 to 100 across genes.

## Filtering rules (raw 31,752 rows -> SOD1 model-ready 1,026 rows)

The following steps are currently only partially documented in code:

1. Selection of `Target_gene == "SOD-1"` rows (2,181 rows -> `training/SOD1_main.csv`): **TODO** — document where/when this filter was applied.
2. De-duplication to one record per (ISIS, Sequence, Location, Smiles), keeping the highest `ASO_volume(nM)` and longest `Treatment_Period(hours)` (2,181 -> 1,026 rows -> `SOD1_cleaned_v1.csv`): implemented in `src/aso_prepare/sod1_pipeline.py::clean_replicate_assays`.
3. Column pruning to produce `SOD1_cleaned_v2.csv` / `SOD1_cleaned_v3.csv` (assay-condition columns dropped): **TODO** — no code currently records this step.
4. Transcript mapping against NCBI RefSeq (fetched FASTA is NM_000454.5): performed in notebook `python/04_group-by_leakage_mapping.ipynb` with a personal NCBI API key; not yet ported to `src/aso_prepare`. **TODO** — document accession version policy and retrieval date.
5. Treatment of implausible inhibition values (e.g. negative percentages present in the raw export for other genes): **TODO** — document any clipping/filtering rule.

## Contact

- Dataset maintainer: **TODO**
