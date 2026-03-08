import os
import sys

import pandas as pd
import pytest

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
scripts_dir = os.path.join(project_root, "workflow/scripts")
for d in [project_root, scripts_dir]:
    if d not in sys.path:
        sys.path.insert(0, d)

from workflow.scripts.classify_orf_confidence import classify_row, main


# ---------------------------------------------------------------------------
# classify_row
# ---------------------------------------------------------------------------

THRESHOLDS = dict(high_min_aa=300, high_min_intrinsic=0.6, putative_min_aa=150, putative_min_intrinsic=0.5)


def test_classify_high_confidence():
    r = {"aa_len": 400, "intrinsic_score": 0.7, "domain_support": True}
    assert classify_row(r, **THRESHOLDS) == "high_confidence_coding"


def test_classify_high_requires_domain_support():
    r = {"aa_len": 400, "intrinsic_score": 0.7, "domain_support": False}
    # meets length+intrinsic but no domain → putative
    assert classify_row(r, **THRESHOLDS) == "putative_coding"


def test_classify_high_requires_intrinsic():
    r = {"aa_len": 400, "intrinsic_score": 0.3, "domain_support": True}
    assert classify_row(r, **THRESHOLDS) == "putative_coding"


def test_classify_putative_by_intrinsic():
    r = {"aa_len": 200, "intrinsic_score": 0.55, "domain_support": False}
    assert classify_row(r, **THRESHOLDS) == "putative_coding"


def test_classify_putative_by_domain():
    r = {"aa_len": 200, "intrinsic_score": 0.1, "domain_support": True}
    assert classify_row(r, **THRESHOLDS) == "putative_coding"


def test_classify_unlikely():
    r = {"aa_len": 100, "intrinsic_score": 0.1, "domain_support": False}
    assert classify_row(r, **THRESHOLDS) == "unlikely_coding"


def test_classify_missing_fields_defaults_to_zero():
    # missing keys default to 0 / False
    assert classify_row({}, **THRESHOLDS) == "unlikely_coding"


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def _write_tsv(path, df):
    df.to_csv(str(path), sep="\t", index=False)


@pytest.fixture
def inputs(tmp_path):
    orfs = pd.DataFrame({
        "te_id": ["te1", "te1", "te2"],
        "orf_id": ["te1|orf_1", "te1|orf_2", "te2|orf_1"],
        "aa_len": [400, 200, 100],
    })
    intrinsic = pd.DataFrame({
        "orf_id": ["te1|orf_1", "te1|orf_2", "te2|orf_1"],
        "intrinsic_score": [0.75, 0.55, 0.1],
        "intrinsic_label": ["high", "medium", "low"],
    })
    domain = pd.DataFrame({
        "orf_id": ["te1|orf_1", "te2|orf_1"],
        "domain_support": [True, False],
        "hit_count": [3, 0],
        "best_domain": ["RT", ""],
        "best_evalue": [1e-10, float("nan")],
        "best_bitscore": [200.0, float("nan")],
    })
    _write_tsv(tmp_path / "orfs.tsv", orfs)
    _write_tsv(tmp_path / "intrinsic.tsv", intrinsic)
    _write_tsv(tmp_path / "domain.tsv", domain)
    return tmp_path


def test_main_output_columns(inputs):
    out = str(inputs / "classified.tsv")
    main(str(inputs / "orfs.tsv"), str(inputs / "intrinsic.tsv"), str(inputs / "domain.tsv"),
         out, **THRESHOLDS)
    df = pd.read_csv(out, sep="\t")
    for col in ("te_id", "orf_id", "aa_len", "intrinsic_score", "domain_support", "confidence_class"):
        assert col in df.columns


def test_main_classification_results(inputs):
    out = str(inputs / "classified.tsv")
    main(str(inputs / "orfs.tsv"), str(inputs / "intrinsic.tsv"), str(inputs / "domain.tsv"),
         out, **THRESHOLDS)
    df = pd.read_csv(out, sep="\t").set_index("orf_id")
    assert df.loc["te1|orf_1", "confidence_class"] == "high_confidence_coding"
    assert df.loc["te1|orf_2", "confidence_class"] == "putative_coding"
    assert df.loc["te2|orf_1", "confidence_class"] == "unlikely_coding"


def test_main_missing_optional_columns_filled(inputs):
    """ORFs with no domain hit (left-join) get NA-filled optional columns."""
    out = str(inputs / "classified.tsv")
    main(str(inputs / "orfs.tsv"), str(inputs / "intrinsic.tsv"), str(inputs / "domain.tsv"),
         out, **THRESHOLDS)
    df = pd.read_csv(out, sep="\t")
    # te1|orf_2 had no domain row → domain_support should be False after fillna
    row = df[df["orf_id"] == "te1|orf_2"].iloc[0]
    assert row["domain_support"] == False  # noqa: E712


def test_main_fills_missing_columns(tmp_path):
    """When domain_df lacks optional columns, they are filled with NA."""
    orfs = pd.DataFrame({"te_id": ["te1"], "orf_id": ["te1|orf_1"], "aa_len": [400]})
    intrinsic = pd.DataFrame({
        "orf_id": ["te1|orf_1"],
        "intrinsic_score": [0.75],
        "intrinsic_label": ["high"],
    })
    # domain_df intentionally missing hit_count / best_domain / best_evalue / best_bitscore
    domain = pd.DataFrame({"orf_id": ["te1|orf_1"], "domain_support": [True]})
    _write_tsv(tmp_path / "orfs.tsv", orfs)
    _write_tsv(tmp_path / "intrinsic.tsv", intrinsic)
    _write_tsv(tmp_path / "domain.tsv", domain)
    out = str(tmp_path / "classified.tsv")
    main(str(tmp_path / "orfs.tsv"), str(tmp_path / "intrinsic.tsv"), str(tmp_path / "domain.tsv"),
         out, **THRESHOLDS)
    df = pd.read_csv(out, sep="\t")
    # optional columns that were absent from domain_df should be present but NA
    assert "best_domain" in df.columns
    assert "hit_count" in df.columns


def test_main_na_intrinsic_score_does_not_raise(tmp_path):
    """ORFs missing intrinsic scores (NA after left join) are classified as unlikely."""
    orfs = pd.DataFrame({"te_id": ["te1"], "orf_id": ["te1|orf_1"], "aa_len": [400]})
    intrinsic = pd.DataFrame(columns=["orf_id", "intrinsic_score", "intrinsic_label"])
    domain = pd.DataFrame({"orf_id": ["te1|orf_1"], "domain_support": [False], "hit_count": [0],
                           "best_domain": [""], "best_evalue": [float("nan")], "best_bitscore": [float("nan")]})
    _write_tsv(tmp_path / "orfs.tsv", orfs)
    _write_tsv(tmp_path / "intrinsic.tsv", intrinsic)
    _write_tsv(tmp_path / "domain.tsv", domain)
    out = str(tmp_path / "classified.tsv")
    main(str(tmp_path / "orfs.tsv"), str(tmp_path / "intrinsic.tsv"), str(tmp_path / "domain.tsv"),
         out, **THRESHOLDS)
    df = pd.read_csv(out, sep="\t")
    assert df.iloc[0]["confidence_class"] == "unlikely_coding"


def test_main_creates_output_dir(tmp_path, inputs):
    out = str(tmp_path / "nested" / "dir" / "classified.tsv")
    main(str(inputs / "orfs.tsv"), str(inputs / "intrinsic.tsv"), str(inputs / "domain.tsv"),
         out, **THRESHOLDS)
    assert os.path.exists(out)
