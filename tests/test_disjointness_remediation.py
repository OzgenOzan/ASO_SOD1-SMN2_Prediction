"""Regression tests for the train/validation disjointness guard (ASO2-1).

These tests use tiny synthetic in-memory inputs only; they never touch the
committed dataset files and never require the optional dependencies.
"""

import inspect

import pandas as pd
import pytest

from src.aso_prepare import sod1_pipeline as sp


def _training_table() -> pd.DataFrame:
    """Two synthetic training rows."""
    return pd.DataFrame(
        {
            "ISIS": ["100001", "100002"],
            "Sequence": ["ACGTACGTACGTACGTACGT", "TTTTCCCCAAAAGGGGTTTT"],
        }
    )


# ---------------------------------------------------------------------------
# validate_no_validation_overlap
# ---------------------------------------------------------------------------

def test_validator_raises_on_overlapping_validation_id():
    df = _training_table()
    with pytest.raises(ValueError, match="ASO2-1"):
        sp.validate_no_validation_overlap(df, validation_ids=["100001"])


def test_validator_raises_on_overlapping_sequence():
    df = _training_table()
    # Sequence comparison is case-insensitive / whitespace-tolerant.
    with pytest.raises(ValueError, match="ASO2-1"):
        sp.validate_no_validation_overlap(df, validation_sequences=[" acgtacgtacgtacgtacgt "])


def test_validator_passes_when_disjoint():
    df = _training_table()
    sp.validate_no_validation_overlap(
        df,
        validation_ids=["999999"],
        validation_sequences=["GGGGGGGGGGGGGGGGGGGG"],
    )


def test_validator_noop_when_no_identifiers_given():
    df = _training_table()
    sp.validate_no_validation_overlap(df)


# ---------------------------------------------------------------------------
# load_validation_identifiers
# ---------------------------------------------------------------------------

def test_load_validation_identifiers_missing_file_warns(tmp_path):
    with pytest.warns(UserWarning, match="Validation CSV not found"):
        ids, sequences = sp.load_validation_identifiers(tmp_path / "absent.csv")
    assert ids is None
    assert sequences is None


def test_load_validation_identifiers_reads_ids_and_sequences(tmp_path):
    csv_path = tmp_path / "validation.csv"
    pd.DataFrame(
        {"ISIS": ["555555"], "Sequence": ["acgtacgtacgtacgtacgt "]}
    ).to_csv(csv_path, index=False)
    ids, sequences = sp.load_validation_identifiers(csv_path)
    assert ids == ["555555"]
    assert sequences == ["ACGTACGTACGTACGTACGT"]


# ---------------------------------------------------------------------------
# Wiring into the export path
# ---------------------------------------------------------------------------

def test_export_path_invokes_validator():
    """Static guard: export_model_ready must call the overlap validator."""
    source = inspect.getsource(sp.export_model_ready)
    assert "validate_no_validation_overlap" in source


def _write_export_inputs(tmp_path, df, validation_rows):
    in_csv = tmp_path / "in.csv"
    df.to_csv(in_csv, index=False)
    validation_csv = tmp_path / "validation.csv"
    pd.DataFrame(validation_rows).to_csv(validation_csv, index=False)
    return in_csv, validation_csv


def test_export_model_ready_raises_on_overlap(tmp_path):
    """End-to-end: export must fail loudly when a validation ID overlaps."""
    df = _training_table()
    df["split_2way"] = ["train", "train"]
    in_csv, validation_csv = _write_export_inputs(
        tmp_path,
        df,
        {"ISIS": ["100002"], "Sequence": ["GGGGGGGGGGGGGGGGGGGG"]},
    )
    with pytest.raises(ValueError, match="ASO2-1"):
        sp.export_model_ready(
            input_csv=in_csv,
            output_dir=tmp_path / "out",
            validation_csv=validation_csv,
        )


def test_export_model_ready_passes_when_disjoint(tmp_path):
    df = _training_table()
    df["split_2way"] = ["train", "test"]
    in_csv, validation_csv = _write_export_inputs(
        tmp_path,
        df,
        {"ISIS": ["999999"], "Sequence": ["GGGGGGGGGGGGGGGGGGGG"]},
    )
    out_dir = tmp_path / "out"
    sp.export_model_ready(input_csv=in_csv, output_dir=out_dir, validation_csv=validation_csv)
    assert (out_dir / "SOD1_model_ready_v1.csv").exists()
