"""Clean, documented SOD1 dataset preparation pipeline.

The original project stores the SOD1 workflow as notebooks. This module keeps the
same scientific assumptions but packages them as rerunnable functions with
project-relative paths.
"""

from __future__ import annotations

import re
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd

from . import paths


SUGAR_MAP = {"M": "MOE", "C": "cEt", "d": "DNA"}
VALID_BASES = set("ACGT")
VALID_PATTERN_SYMBOLS = set(SUGAR_MAP)


def clean_replicate_assays(
    input_csv: Path = paths.SOD1_RAW,
    output_csv: Path = paths.SOD1_CLEANED_V1,
) -> pd.DataFrame:
    """Keep the highest-dose, longest-duration row per ASO target record.

    The raw export can contain repeated rows for the same `ISIS`, `Sequence`,
    `Location`, and `Smiles`. The original notebook sorted dose and duration in
    descending order, then kept the first row per group. This function preserves
    that rule and writes `SOD1_cleaned_v1.csv`.
    """

    df = pd.read_csv(input_csv)
    df.columns = df.columns.str.strip()

    sort_columns = [
        "ISIS",
        "Sequence",
        "Location",
        "Smiles",
        "ASO_volume(nM)",
        "Treatment_Period(hours)",
    ]
    dedupe_columns = ["ISIS", "Sequence", "Location", "Smiles"]

    missing = sorted(set(sort_columns) - set(df.columns))
    if missing:
        raise ValueError(f"Missing columns required for assay cleanup: {missing}")

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    cleaned = (
        df.sort_values(
            by=sort_columns,
            ascending=[True, True, True, True, False, False],
        )
        .drop_duplicates(subset=dedupe_columns, keep="first")
        .reset_index(drop=True)
    )
    cleaned.to_csv(output_csv, index=False)
    return cleaned


def validate_lengths(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Normalize ASO sequence fields and flag sequence/chemistry length issues."""

    df = df.copy()
    df["Sequence"] = df["Sequence"].astype(str).str.upper().str.strip()
    df["Chemical_Pattern"] = df["Chemical_Pattern"].astype(str).str.strip()

    if "seq_length" in df.columns and "seq_length_original" not in df.columns:
        df["seq_length_original"] = df["seq_length"]

    df["seq_length_recalc"] = df["Sequence"].str.len()
    df["chem_pattern_length"] = df["Chemical_Pattern"].str.len()
    df["seq_length"] = df["seq_length_recalc"]
    df["valid_sequence"] = df["Sequence"].apply(lambda value: set(value).issubset(VALID_BASES))
    df["valid_chemical_pattern"] = df["Chemical_Pattern"].apply(
        lambda value: set(value).issubset(VALID_PATTERN_SYMBOLS)
    )
    df["lengths_match"] = (
        (df["seq_length_recalc"] == df["chem_pattern_length"])
        & df["valid_sequence"]
        & df["valid_chemical_pattern"]
    )

    issue_columns = [
        "row_id",
        "ISIS",
        "Sequence",
        "Chemical_Pattern",
        "seq_length_original",
        "seq_length_recalc",
        "chem_pattern_length",
        "valid_sequence",
        "valid_chemical_pattern",
        "lengths_match",
    ]
    issue_columns = [col for col in issue_columns if col in df.columns]
    issues = df.loc[~df["lengths_match"], issue_columns]
    return df, issues


def build_position_annotations(df: pd.DataFrame) -> pd.DataFrame:
    """Expand each ASO into one row per nucleotide and annotate sugar chemistry."""

    records: list[dict[str, object]] = []
    for _, row in df.iterrows():
        sequence = row["Sequence"]
        pattern = row["Chemical_Pattern"]
        if len(sequence) != len(pattern):
            raise ValueError(f"Length mismatch for ISIS={row['ISIS']}")

        for index, (base, symbol) in enumerate(zip(sequence, pattern)):
            sugar = SUGAR_MAP.get(symbol, "UNKNOWN")
            records.append(
                {
                    "row_id": row["row_id"],
                    "ISIS": row["ISIS"],
                    "position_0based": index,
                    "position_1based": index + 1,
                    "base": base,
                    "chem_symbol": symbol,
                    "sugar": sugar,
                    "is_DNA_gap": sugar == "DNA",
                    "is_MOE": sugar == "MOE",
                    "is_cEt": sugar == "cEt",
                    "is_modified_sugar": sugar in {"MOE", "cEt"},
                    "seq_length": len(sequence),
                }
            )
    return pd.DataFrame(records)


def parse_linkage_locations(location: object) -> set[int]:
    """Extract 1-based linkage positions from strings such as `1?2?3/else`."""

    if pd.isna(location):
        return set()
    text = str(location).strip()
    if text == "" or text.lower() == "nan" or text == "else":
        return set()
    first_part = text.split("/")[0].strip()
    if first_part == "" or first_part == "else":
        return set()
    return {int(value) for value in re.findall(r"\d+", first_part)}


def build_linkage_annotations(df: pd.DataFrame) -> pd.DataFrame:
    """Expand each ASO into one row per internucleotide linkage.

    Assumptions inherited from the notebooks:
    `phosphorothioate` means all linkages are PS, while
    `phosphodiester/phosphorothioate` means listed positions are PO and the rest
    are PS. Listed positions are interpreted as 1-based.
    """

    records: list[dict[str, object]] = []
    for _, row in df.iterrows():
        listed_positions = parse_linkage_locations(row["Linkage_Location"])
        for linkage_position in range(1, int(row["seq_length"])):
            linkage_type = str(row["Linkage"]).strip()
            if linkage_type == "phosphorothioate":
                linkage = "PS"
            elif linkage_type == "phosphodiester/phosphorothioate":
                linkage = "PO" if linkage_position in listed_positions else "PS"
            else:
                linkage = "UNKNOWN"

            records.append(
                {
                    "row_id": row["row_id"],
                    "ISIS": row["ISIS"],
                    "linkage_position_1based": linkage_position,
                    "between_nt_1based": f"{linkage_position}-{linkage_position + 1}",
                    "linkage": linkage,
                    "is_PS": linkage == "PS",
                    "is_PO": linkage == "PO",
                    "is_unknown_linkage": linkage == "UNKNOWN",
                }
            )
    return pd.DataFrame(records)


def add_chemistry_summaries(
    df: pd.DataFrame,
    position_df: pd.DataFrame,
    linkage_df: pd.DataFrame,
) -> pd.DataFrame:
    """Add per-ASO counts and fractions derived from chemistry annotations."""

    linkage_summary = (
        linkage_df.groupby("row_id")
        .agg(
            n_linkages=("linkage", "size"),
            n_PS=("is_PS", "sum"),
            n_PO=("is_PO", "sum"),
            n_unknown_linkage=("is_unknown_linkage", "sum"),
        )
        .reset_index()
    )
    chem_summary = (
        position_df.groupby("row_id")
        .agg(
            n_MOE=("is_MOE", "sum"),
            n_cEt=("is_cEt", "sum"),
            n_DNA=("is_DNA_gap", "sum"),
            n_modified_sugar=("is_modified_sugar", "sum"),
        )
        .reset_index()
    )

    df = df.merge(linkage_summary, on="row_id", how="left", validate="one_to_one")
    df = df.merge(chem_summary, on="row_id", how="left", validate="one_to_one")
    df["expected_n_linkages"] = df["seq_length"] - 1
    df["linkage_count_valid"] = df["n_linkages"] == df["expected_n_linkages"]
    df["frac_MOE"] = df["n_MOE"] / df["seq_length"]
    df["frac_cEt"] = df["n_cEt"] / df["seq_length"]
    df["frac_DNA"] = df["n_DNA"] / df["seq_length"]
    df["frac_modified_sugar"] = df["n_modified_sugar"] / df["seq_length"]
    df["chem_count_sum_valid"] = (df["n_MOE"] + df["n_cEt"] + df["n_DNA"]) == df["seq_length"]
    return df


def build_chemistry_stage(
    input_csv: Path = paths.SOD1_CLEANED_V3,
    output_dir: Path = paths.SOD1_CHEMISTRY_DIR,
) -> pd.DataFrame:
    """Create SOD1 chemistry, position-level, and linkage-level CSV outputs."""

    output_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(input_csv)
    df.columns = df.columns.str.strip()

    required = ["ISIS", "Sequence", "Chemical_Pattern", "Linkage", "Linkage_Location"]
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise ValueError(f"Missing columns required for chemistry stage: {missing}")

    df = df.reset_index(drop=True)
    df["row_id"] = np.arange(len(df))
    df, length_issues = validate_lengths(df)
    if len(length_issues) > 0:
        length_issues.to_csv(output_dir / "SOD1_validation_issues_v1.csv", index=False)
        raise ValueError("Length or symbol validation failed.")

    position_df = build_position_annotations(df)
    linkage_df = build_linkage_annotations(df)
    df = add_chemistry_summaries(df, position_df, linkage_df)

    checks = ["lengths_match", "linkage_count_valid", "chem_count_sum_valid"]
    if not df[checks].all().all() or df["n_unknown_linkage"].sum() > 0:
        raise ValueError("Chemistry validation failed.")

    df.to_csv(output_dir / "SOD1_process_v1.csv", index=False)
    position_df.to_csv(output_dir / "SOD1_position_level_annotations_v1.csv", index=False)
    linkage_df.to_csv(output_dir / "SOD1_linkage_level_annotations_v1.csv", index=False)
    write_data_dictionary(df, output_dir / "SOD1_process_v1_data_dictionary.csv")
    return df


def add_inhibition_groups_and_splits(
    input_csv: Path = paths.SOD1_CHEMISTRY_DIR / "SOD1_process_v1.csv",
    output_dir: Path = paths.SOD1_SPLIT_DIR,
) -> pd.DataFrame:
    """Create inhibition quartile labels and sequence-grouped train/test splits."""

    from sklearn.model_selection import StratifiedGroupKFold

    output_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(input_csv)
    df.columns = df.columns.str.strip()
    df["Sequence"] = df["Sequence"].astype(str).str.upper().str.strip()
    df["inhibition_pct"] = pd.to_numeric(df["Inhibition(%)"], errors="coerce")
    df["inhibition_rank_first"] = df["inhibition_pct"].rank(method="first")
    df["inhibition_group_4_equal_count"] = pd.qcut(
        df["inhibition_rank_first"],
        q=4,
        labels=["low", "mid_low", "mid_high", "high"],
    )

    sgkf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)
    train_index, test_index = next(
        sgkf.split(
            np.zeros(len(df)),
            df["inhibition_group_4_equal_count"].astype(str),
            groups=df["Sequence"],
        )
    )
    df["split_2way"] = "train"
    df.loc[test_index, "split_2way"] = "test"

    group_summary = (
        df.groupby("inhibition_group_4_equal_count", observed=False)
        .agg(
            n_ASOs=("ISIS", "count"),
            min_inhibition=("inhibition_pct", "min"),
            max_inhibition=("inhibition_pct", "max"),
            mean_inhibition=("inhibition_pct", "mean"),
            median_inhibition=("inhibition_pct", "median"),
        )
        .reset_index()
    )
    df.to_csv(output_dir / "SOD1_process_v2_inhibition_groups_and_splits.csv", index=False)
    group_summary.to_csv(output_dir / "SOD1_inhibition_group_summary_v2.csv", index=False)
    return df


def add_sequence_features(df: pd.DataFrame) -> pd.DataFrame:
    """Add base composition and simple run-length features."""

    df = df.copy()
    sequence = df["Sequence"].astype(str).str.upper().str.strip()
    df["gc_fraction"] = sequence.apply(lambda value: (value.count("G") + value.count("C")) / len(value))
    for base in "ACGT":
        df[f"n_{base}"] = sequence.str.count(base)
        df[f"frac_{base}"] = df[f"n_{base}"] / df["seq_length"]
    df["longest_homopolymer"] = sequence.apply(
        lambda value: max((len(match.group(0)) for match in re.finditer(r"(A+|C+|G+|T+)", value)), default=0)
    )
    df["has_homopolymer_4plus"] = df["longest_homopolymer"] >= 4
    df["has_homopolymer_5plus"] = df["longest_homopolymer"] >= 5
    df["max_gc_stretch"] = sequence.apply(
        lambda value: max((len(match.group(0)) for match in re.finditer(r"[GC]+", value)), default=0)
    )
    df["base_count_sum_valid"] = df[[f"n_{base}" for base in "ACGT"]].sum(axis=1) == df["seq_length"]
    return df


def add_gap_features(df: pd.DataFrame) -> pd.DataFrame:
    """Describe the contiguous DNA gap encoded by `d` in `Chemical_Pattern`."""

    def dna_gap_bounds(pattern: str) -> tuple[int, int, int]:
        matches = list(re.finditer(r"d+", str(pattern)))
        if not matches:
            return -1, -1, 0
        longest = max(matches, key=lambda match: len(match.group(0)))
        return longest.start(), longest.end(), len(longest.group(0))

    df = df.copy()
    gap_info = df["Chemical_Pattern"].apply(dna_gap_bounds).apply(pd.Series)
    gap_info.columns = ["gap_start_0based", "gap_end_0based_exclusive", "gap_length"]
    df = pd.concat([df.reset_index(drop=True), gap_info.reset_index(drop=True)], axis=1)
    df["longest_DNA_gap"] = df["gap_length"]
    df["gap_sequence"] = df.apply(
        lambda row: row["Sequence"][row["gap_start_0based"] : row["gap_end_0based_exclusive"]]
        if row["gap_length"] > 0
        else "",
        axis=1,
    )
    df["flank5_length"] = df["gap_start_0based"].clip(lower=0)
    df["flank3_length"] = df["seq_length"] - df["gap_end_0based_exclusive"].clip(lower=0)
    df["flank_asymmetry_5minus3"] = df["flank5_length"] - df["flank3_length"]
    df["gap_length_valid"] = df["gap_length"] == df["Chemical_Pattern"].str.count("d")
    return df


def add_gap_3mer_features(df: pd.DataFrame) -> pd.DataFrame:
    """Count every possible DNA 3-mer inside the gap sequence."""

    kmers = ["".join(kmer) for kmer in product("ACGT", repeat=3)]

    def count_kmers(gap_sequence: object) -> dict[str, int]:
        text = str(gap_sequence).upper()
        counts = {f"gap_3mer_{kmer}": 0 for kmer in kmers}
        for index in range(max(len(text) - 2, 0)):
            key = f"gap_3mer_{text[index:index + 3]}"
            if key in counts:
                counts[key] += 1
        return counts

    kmer_df = pd.DataFrame(df["gap_sequence"].apply(count_kmers).tolist())
    df = pd.concat([df.reset_index(drop=True), kmer_df.reset_index(drop=True)], axis=1)
    df["gap_3mer_total_count"] = kmer_df.sum(axis=1)
    df["gap_3mer_expected_count"] = df["gap_length"].apply(lambda value: max(int(value) - 2, 0))
    df["gap_3mer_count_valid"] = df["gap_3mer_total_count"] == df["gap_3mer_expected_count"]
    return df


def add_transcript_coordinate_features(
    df: pd.DataFrame,
    fasta_path: Path = paths.SOD1_TRANSCRIPT_DIR / "SOD1_refseq_transcript_fetched_from_NCBI.fasta",
) -> pd.DataFrame:
    """Add target-position features derived from transcript mapping columns.

    These features summarize where the ASO binds within the fetched SOD1
    transcript. Context-window sequences are added when the RefSeq FASTA file is
    available and Biopython can parse it.
    """

    df = df.copy()
    required = [
        "n_transcript_hits",
        "target_start_1based",
        "target_end_1based_inclusive",
        "transcript_length",
    ]
    if sorted(set(required) - set(df.columns)):
        df["target_midpoint_1based"] = pd.NA
        df["target_relative_start"] = pd.NA
        df["target_relative_end"] = pd.NA
        df["target_relative_midpoint"] = pd.NA
        df["target_has_exact_refseq_match"] = False
        df["target_mapping_unique"] = False
        df["target_mapping_multiple"] = False
    else:
        n_hits = pd.to_numeric(df["n_transcript_hits"], errors="coerce").fillna(0)
        target_start = pd.to_numeric(df["target_start_1based"], errors="coerce")
        target_end = pd.to_numeric(df["target_end_1based_inclusive"], errors="coerce")
        transcript_length = pd.to_numeric(df["transcript_length"], errors="coerce")
        midpoint = (target_start + target_end) / 2

        df["target_midpoint_1based"] = midpoint
        df["target_relative_start"] = target_start / transcript_length
        df["target_relative_end"] = target_end / transcript_length
        df["target_relative_midpoint"] = midpoint / transcript_length
        df["target_has_exact_refseq_match"] = n_hits.astype(int) > 0
        df["target_mapping_unique"] = n_hits.astype(int) == 1
        df["target_mapping_multiple"] = n_hits.astype(int) > 1

    for flank in [5, 10, 20]:
        df[f"target_context_{flank}nt_each_side"] = pd.NA

    if not fasta_path.exists():
        return df

    try:
        from Bio import SeqIO
    except ImportError:
        return df

    record = next(SeqIO.parse(str(fasta_path), "fasta"), None)
    if record is None:
        return df

    transcript_sequence = str(record.seq).upper().replace("U", "T")

    def extract_context(row: pd.Series, flank: int) -> object:
        start = row.get("target_start_1based", pd.NA)
        end = row.get("target_end_1based_inclusive", pd.NA)
        if pd.isna(start) or pd.isna(end):
            return pd.NA
        start_0based = int(start) - 1
        end_exclusive = int(end)
        context_start = max(0, start_0based - flank)
        context_end = min(len(transcript_sequence), end_exclusive + flank)
        return transcript_sequence[context_start:context_end]

    for flank in [5, 10, 20]:
        df[f"target_context_{flank}nt_each_side"] = df.apply(
            lambda row: extract_context(row, flank),
            axis=1,
        )
    return df


def add_basic_tm_features(df: pd.DataFrame) -> pd.DataFrame:
    """Add optional Biopython nearest-neighbor DNA Tm estimates."""

    df = df.copy()
    try:
        from Bio.SeqUtils import MeltingTemp as mt
    except ImportError:
        df["basic_DNA_Tm_NN_full_sequence"] = np.nan
        df["basic_DNA_Tm_NN_gap_sequence"] = np.nan
        return df

    def tm_nn(sequence: object) -> float:
        text = str(sequence).upper()
        if len(text) < 2:
            return np.nan
        try:
            return float(mt.Tm_NN(text, nn_table=mt.DNA_NN4))
        except Exception:
            return np.nan

    df["basic_DNA_Tm_NN_full_sequence"] = df["Sequence"].apply(tm_nn)
    df["basic_DNA_Tm_NN_gap_sequence"] = df["gap_sequence"].apply(tm_nn)
    return df


def add_rdkit_features(df: pd.DataFrame) -> pd.DataFrame:
    """Add optional molecular descriptors from the `Smiles` column using RDKit."""

    rdkit_columns = [
        "canonical_smiles",
        "rdkit_valid_smiles",
        "rdkit_mol_wt",
        "rdkit_heavy_atom_count",
        "rdkit_num_atoms",
        "rdkit_num_bonds",
        "rdkit_num_rings",
        "rdkit_tpsa",
        "rdkit_num_h_donors",
        "rdkit_num_h_acceptors",
        "rdkit_formal_charge",
    ]
    df = df.copy()
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors
    except ImportError:
        for column in rdkit_columns:
            df[column] = pd.NA
        return df

    def descriptors(smiles: object) -> dict[str, object]:
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            return {column: np.nan for column in rdkit_columns} | {"rdkit_valid_smiles": False}
        return {
            "canonical_smiles": Chem.MolToSmiles(mol, canonical=True),
            "rdkit_valid_smiles": True,
            "rdkit_mol_wt": Descriptors.MolWt(mol),
            "rdkit_heavy_atom_count": Descriptors.HeavyAtomCount(mol),
            "rdkit_num_atoms": mol.GetNumAtoms(),
            "rdkit_num_bonds": mol.GetNumBonds(),
            "rdkit_num_rings": rdMolDescriptors.CalcNumRings(mol),
            "rdkit_tpsa": rdMolDescriptors.CalcTPSA(mol),
            "rdkit_num_h_donors": rdMolDescriptors.CalcNumHBD(mol),
            "rdkit_num_h_acceptors": rdMolDescriptors.CalcNumHBA(mol),
            "rdkit_formal_charge": Chem.GetFormalCharge(mol),
        }

    descriptor_df = pd.DataFrame(df["Smiles"].apply(descriptors).tolist())
    return pd.concat([df.reset_index(drop=True), descriptor_df.reset_index(drop=True)], axis=1)


def build_feature_stage(
    input_csv: Path = paths.SOD1_TRANSCRIPT_DIR / "SOD1_process_v3_transcript_mapping.csv",
    output_dir: Path = paths.SOD1_FEATURE_DIR,
) -> pd.DataFrame:
    """Build the non-modeling feature table from mapped SOD1 rows."""

    output_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(input_csv)
    df.columns = df.columns.str.strip()
    df = add_sequence_features(df)
    df.to_csv(output_dir / "SOD1_process_v4a_sequence_features.csv", index=False)
    df = add_gap_features(df)
    df.to_csv(output_dir / "SOD1_process_v4b_gap_architecture_features.csv", index=False)
    df = add_gap_3mer_features(df)
    df.to_csv(output_dir / "SOD1_process_v4c_gap_3mer_features.csv", index=False)
    df = add_transcript_coordinate_features(df)
    df.to_csv(output_dir / "SOD1_process_v4d_transcript_coordinate_features.csv", index=False)
    df = add_basic_tm_features(df)
    df.to_csv(output_dir / "SOD1_process_v4e_basic_tm_features.csv", index=False)
    df = add_rdkit_features(df)
    df.to_csv(output_dir / "SOD1_process_v4f_rdkit_features.csv", index=False)
    return df


def write_data_dictionary(df: pd.DataFrame, output_csv: Path) -> None:
    """Write a compact schema summary for a generated dataset."""

    dictionary = pd.DataFrame(
        {
            "column": df.columns,
            "dtype": [str(df[column].dtype) for column in df.columns],
            "n_missing": [df[column].isna().sum() for column in df.columns],
            "n_unique": [df[column].nunique(dropna=True) for column in df.columns],
        }
    )
    dictionary.to_csv(output_csv, index=False)


def export_model_ready(
    input_csv: Path = paths.SOD1_FEATURE_DIR / "SOD1_process_v4f_rdkit_features.csv",
    output_dir: Path = paths.SOD1_MODEL_READY_DIR,
) -> pd.DataFrame:
    """Validate and export the final model-ready SOD1 feature table."""

    output_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(input_csv)
    df.columns = df.columns.str.strip()

    checks: dict[str, bool] = {}
    for column in [
        "lengths_match",
        "valid_sequence",
        "valid_chemical_pattern",
        "linkage_count_valid",
        "chem_count_sum_valid",
        "base_count_sum_valid",
        "gap_length_valid",
        "gap_3mer_count_valid",
    ]:
        if column in df.columns:
            checks[f"{column}_all_true"] = bool(df[column].all())

    if "split_2way" in df.columns:
        train_sequences = set(df.loc[df["split_2way"] == "train", "Sequence"])
        test_sequences = set(df.loc[df["split_2way"] == "test", "Sequence"])
        checks["no_sequence_leakage_train_test_2way"] = not train_sequences.intersection(test_sequences)

    validation_summary = pd.DataFrame({"check": list(checks), "passed": list(checks.values())})
    df.to_csv(output_dir / "SOD1_model_ready_v1.csv", index=False)
    write_data_dictionary(df, output_dir / "SOD1_model_ready_v1_data_dictionary.csv")
    validation_summary.to_csv(output_dir / "SOD1_model_ready_v1_validation_summary.csv", index=False)
    return df


def run_sod1_pipeline() -> pd.DataFrame:
    """Run the documented SOD1 pipeline from chemistry annotations to final export."""

    build_chemistry_stage()
    add_inhibition_groups_and_splits()
    if (paths.SOD1_TRANSCRIPT_DIR / "SOD1_process_v3_transcript_mapping.csv").exists():
        build_feature_stage()
        return export_model_ready()
    return pd.read_csv(paths.SOD1_SPLIT_DIR / "SOD1_process_v2_inhibition_groups_and_splits.csv")


if __name__ == "__main__":
    result = run_sod1_pipeline()
    print(f"SOD1 pipeline complete: {result.shape[0]} rows, {result.shape[1]} columns")
