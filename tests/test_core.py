"""Smoke tests for the core SOD1 preparation helpers in src/aso_prepare.

These tests use tiny synthetic inputs only; they never touch the committed
dataset files and never require the optional dependencies (Biopython/RDKit).
"""

import numpy as np
import pandas as pd
import pytest

from src.aso_prepare import sod1_pipeline as sp


# ---------------------------------------------------------------------------
# parse_linkage_locations
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        (None, set()),
        (np.nan, set()),
        ("", set()),
        ("nan", set()),
        ("else", set()),
        ("/else", set()),
        ("1?2?3/else", {1, 2, 3}),
        ("5", {5}),
        ("10?11?12/C/else", {10, 11, 12}),
    ],
)
def test_parse_linkage_locations(raw, expected):
    assert sp.parse_linkage_locations(raw) == expected


# ---------------------------------------------------------------------------
# Gap architecture (dna_gap_bounds logic via add_gap_features)
# ---------------------------------------------------------------------------

def _gap_df(patterns, sequences):
    return pd.DataFrame(
        {
            "Chemical_Pattern": patterns,
            "Sequence": sequences,
            "seq_length": [len(s) for s in sequences],
        }
    )


def test_gap_features_single_contiguous_gap():
    df = sp.add_gap_features(_gap_df(["MMMMMddddddddddMMMMM"], ["AAAACCCCCCCCCCCAAAA"]))
    row = df.iloc[0]
    assert row["gap_start_0based"] == 5
    assert row["gap_end_0based_exclusive"] == 15
    assert row["gap_length"] == 10
    assert row["gap_sequence"] == "CCCCCCCCCC"
    assert row["flank5_length"] == 5
    assert row["flank3_length"] == 4
    assert bool(row["gap_length_valid"])


def test_gap_features_no_gap():
    df = sp.add_gap_features(_gap_df(["MMMM"], ["ACGT"]))
    row = df.iloc[0]
    assert row["gap_start_0based"] == -1
    assert row["gap_length"] == 0
    assert row["gap_sequence"] == ""
    assert row["flank5_length"] == 0
    assert bool(row["gap_length_valid"])


def test_gap_features_picks_longest_gap_block():
    # Two DNA blocks: longest (3 d's) must win.
    df = sp.add_gap_features(_gap_df(["MddMMdddM"], ["ACGTACGTA"]))
    row = df.iloc[0]
    assert row["gap_start_0based"] == 5
    assert row["gap_length"] == 3
    assert row["gap_sequence"] == "CGT"
    # gap_length_valid compares against the TOTAL d count, so a split gap
    # is expected to fail validation - lock in that documented behaviour.
    assert not bool(row["gap_length_valid"])


# ---------------------------------------------------------------------------
# 3-mer counting (via add_gap_3mer_features)
# ---------------------------------------------------------------------------

def test_gap_3mer_counts_all_unique():
    df = pd.DataFrame({"gap_sequence": ["ACGT"], "gap_length": [4]})
    out = sp.add_gap_3mer_features(df)
    row = out.iloc[0]
    assert row["gap_3mer_ACG"] == 1
    assert row["gap_3mer_CGT"] == 1
    assert row["gap_3mer_total_count"] == 2
    assert row["gap_3mer_expected_count"] == 2
    assert bool(row["gap_3mer_count_valid"])


def test_gap_3mer_counts_overlapping_repeats():
    df = pd.DataFrame({"gap_sequence": ["AAAAA"], "gap_length": [5]})
    out = sp.add_gap_3mer_features(df)
    assert out.iloc[0]["gap_3mer_AAA"] == 3
    assert bool(out.iloc[0]["gap_3mer_count_valid"])


def test_gap_3mer_short_gap_counts_nothing():
    df = pd.DataFrame({"gap_sequence": ["AC", ""], "gap_length": [2, 0]})
    out = sp.add_gap_3mer_features(df)
    assert out["gap_3mer_total_count"].tolist() == [0, 0]
    assert out["gap_3mer_count_valid"].all()


# ---------------------------------------------------------------------------
# Train/test sequence disjointness check (via export_model_ready)
# ---------------------------------------------------------------------------

_TRUE_FLAG_COLUMNS = [
    "lengths_match",
    "valid_sequence",
    "valid_chemical_pattern",
    "linkage_count_valid",
    "chem_count_sum_valid",
    "base_count_sum_valid",
    "gap_length_valid",
    "gap_3mer_count_valid",
]


def _model_ready_input(sequences, splits):
    df = pd.DataFrame({"Sequence": sequences, "split_2way": splits})
    for col in _TRUE_FLAG_COLUMNS:
        df[col] = True
    return df


def _run_export(df, tmp_path):
    in_csv = tmp_path / "in.csv"
    out_dir = tmp_path / "out"
    df.to_csv(in_csv, index=False)
    sp.export_model_ready(input_csv=in_csv, output_dir=out_dir)
    summary = pd.read_csv(out_dir / "SOD1_model_ready_v1_validation_summary.csv")
    return dict(zip(summary["check"], summary["passed"]))


def test_split_disjoint_sequences_passes(tmp_path):
    df = _model_ready_input(["ACGT", "ACGA", "TTTT"], ["train", "train", "test"])
    checks = _run_export(df, tmp_path)
    assert checks["no_sequence_leakage_train_test_2way"] == True  # noqa: E712
    assert all(checks[f"{c}_all_true"] == True for c in _TRUE_FLAG_COLUMNS)  # noqa: E712


def test_split_shared_sequence_is_flagged(tmp_path):
    df = _model_ready_input(["ACGT", "ACGT"], ["train", "test"])
    checks = _run_export(df, tmp_path)
    assert checks["no_sequence_leakage_train_test_2way"] == False  # noqa: E712


def test_committed_model_ready_split_is_sequence_disjoint():
    """Regression guard: the shipped dataset must keep train/test sequences disjoint."""
    from src.aso_prepare import paths

    csv_path = paths.SOD1_MODEL_READY_DIR / "SOD1_model_ready_v1.csv"
    if not csv_path.exists():
        pytest.skip("committed model-ready dataset not present")
    df = pd.read_csv(csv_path)
    train = set(df.loc[df["split_2way"] == "train", "Sequence"])
    test = set(df.loc[df["split_2way"] == "test", "Sequence"])
    assert not train.intersection(test)
