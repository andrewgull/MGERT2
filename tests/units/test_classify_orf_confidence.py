import os
import sys

import pandas as pd
import pytest

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
scripts_dir = os.path.join(project_root, "workflow/scripts")
for d in [project_root, scripts_dir]:
    if d not in sys.path:
        sys.path.insert(0, d)

from workflow.scripts.classify_orf_confidence import (
    apply_domain_filter,
    classify_row,
    main,
)


# ---------------------------------------------------------------------------
# classify_row
# ---------------------------------------------------------------------------

THRESHOLDS = dict(
    high_min_aa=300,
    high_min_intrinsic=0.6,
    putative_min_aa=150,
    putative_min_intrinsic=0.5,
)


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
    orfs = pd.DataFrame(
        {
            "te_id": ["te1", "te1", "te2"],
            "orf_id": ["te1|orf_1", "te1|orf_2", "te2|orf_1"],
            "aa_len": [400, 200, 100],
        }
    )
    intrinsic = pd.DataFrame(
        {
            "orf_id": ["te1|orf_1", "te1|orf_2", "te2|orf_1"],
            "intrinsic_score": [0.75, 0.55, 0.1],
            "intrinsic_label": ["high", "medium", "low"],
        }
    )
    domain = pd.DataFrame(
        {
            "orf_id": ["te1|orf_1", "te2|orf_1"],
            "domain_support": [True, False],
            "hit_count": [3, 0],
            "best_domain": ["RT", ""],
            "best_evalue": [1e-10, float("nan")],
            "best_bitscore": [200.0, float("nan")],
            "domain_types": ["RT", ""],
            "n_domain_types": [1, 0],
        }
    )
    _write_tsv(tmp_path / "orfs.tsv", orfs)
    _write_tsv(tmp_path / "intrinsic.tsv", intrinsic)
    _write_tsv(tmp_path / "domain.tsv", domain)
    return tmp_path


def test_main_output_columns(inputs):
    out = str(inputs / "classified.tsv")
    main(
        str(inputs / "orfs.tsv"),
        str(inputs / "intrinsic.tsv"),
        str(inputs / "domain.tsv"),
        out,
        **THRESHOLDS
    )
    df = pd.read_csv(out, sep="\t")
    for col in (
        "te_id",
        "orf_id",
        "aa_len",
        "intrinsic_score",
        "domain_support",
        "confidence_class",
        "domain_types",
        "n_domain_types",
        "passes_domain_filter",
        "te_domain_complete",
    ):
        assert col in df.columns


def test_main_classification_results(inputs):
    out = str(inputs / "classified.tsv")
    main(
        str(inputs / "orfs.tsv"),
        str(inputs / "intrinsic.tsv"),
        str(inputs / "domain.tsv"),
        out,
        **THRESHOLDS
    )
    df = pd.read_csv(out, sep="\t").set_index("orf_id")
    assert df.loc["te1|orf_1", "confidence_class"] == "high_confidence_coding"
    assert df.loc["te1|orf_2", "confidence_class"] == "putative_coding"
    assert df.loc["te2|orf_1", "confidence_class"] == "unlikely_coding"


def test_main_missing_optional_columns_filled(inputs):
    """ORFs with no domain hit (left-join) get NA-filled optional columns."""
    out = str(inputs / "classified.tsv")
    main(
        str(inputs / "orfs.tsv"),
        str(inputs / "intrinsic.tsv"),
        str(inputs / "domain.tsv"),
        out,
        **THRESHOLDS
    )
    df = pd.read_csv(out, sep="\t")
    # te1|orf_2 had no domain row → domain_support should be False after fillna
    row = df[df["orf_id"] == "te1|orf_2"].iloc[0]
    assert row["domain_support"] == False  # noqa: E712


def test_main_fills_missing_columns(tmp_path):
    """When domain_df lacks optional columns, they are filled with NA."""
    orfs = pd.DataFrame({"te_id": ["te1"], "orf_id": ["te1|orf_1"], "aa_len": [400]})
    intrinsic = pd.DataFrame(
        {
            "orf_id": ["te1|orf_1"],
            "intrinsic_score": [0.75],
            "intrinsic_label": ["high"],
        }
    )
    # domain_df intentionally missing hit_count / best_domain / best_evalue / best_bitscore
    domain = pd.DataFrame({"orf_id": ["te1|orf_1"], "domain_support": [True]})
    _write_tsv(tmp_path / "orfs.tsv", orfs)
    _write_tsv(tmp_path / "intrinsic.tsv", intrinsic)
    _write_tsv(tmp_path / "domain.tsv", domain)
    out = str(tmp_path / "classified.tsv")
    main(
        str(tmp_path / "orfs.tsv"),
        str(tmp_path / "intrinsic.tsv"),
        str(tmp_path / "domain.tsv"),
        out,
        **THRESHOLDS
    )
    df = pd.read_csv(out, sep="\t")
    # optional columns that were absent from domain_df should be present but NA
    assert "best_domain" in df.columns
    assert "hit_count" in df.columns


def test_main_na_intrinsic_score_does_not_raise(tmp_path):
    """ORFs missing intrinsic scores (NA after left join) are classified as unlikely."""
    orfs = pd.DataFrame({"te_id": ["te1"], "orf_id": ["te1|orf_1"], "aa_len": [400]})
    intrinsic = pd.DataFrame(columns=["orf_id", "intrinsic_score", "intrinsic_label"])
    domain = pd.DataFrame(
        {
            "orf_id": ["te1|orf_1"],
            "domain_support": [False],
            "hit_count": [0],
            "best_domain": [""],
            "best_evalue": [float("nan")],
            "best_bitscore": [float("nan")],
        }
    )
    _write_tsv(tmp_path / "orfs.tsv", orfs)
    _write_tsv(tmp_path / "intrinsic.tsv", intrinsic)
    _write_tsv(tmp_path / "domain.tsv", domain)
    out = str(tmp_path / "classified.tsv")
    main(
        str(tmp_path / "orfs.tsv"),
        str(tmp_path / "intrinsic.tsv"),
        str(tmp_path / "domain.tsv"),
        out,
        **THRESHOLDS
    )
    df = pd.read_csv(out, sep="\t")
    assert df.iloc[0]["confidence_class"] == "unlikely_coding"


def test_main_creates_output_dir(tmp_path, inputs):
    out = str(tmp_path / "nested" / "dir" / "classified.tsv")
    main(
        str(inputs / "orfs.tsv"),
        str(inputs / "intrinsic.tsv"),
        str(inputs / "domain.tsv"),
        out,
        **THRESHOLDS
    )
    assert os.path.exists(out)


# ---------------------------------------------------------------------------
# apply_domain_filter
# ---------------------------------------------------------------------------


def _make_df(rows):
    """Build a minimal merged DataFrame for apply_domain_filter tests."""
    return pd.DataFrame(rows)


def test_no_required_domains_te_domain_complete_always_true():
    df = _make_df(
        [
            {"te_id": "te1", "orf_id": "te1|orf_1", "domain_types": ""},
            {"te_id": "te1", "orf_id": "te1|orf_2", "domain_types": "RT"},
        ]
    )
    result = apply_domain_filter(
        df,
        required_domains=[],
        min_domain_types_per_orf=1,
        allow_split_across_orfs=False,
        te_coverage_scope="all",
    )
    assert result["te_domain_complete"].all()


def test_per_orf_filter_single_required_domain():
    df = _make_df(
        [
            {"te_id": "te1", "orf_id": "te1|orf_1", "domain_types": "RT"},
            {"te_id": "te1", "orf_id": "te1|orf_2", "domain_types": "EN"},
            {"te_id": "te2", "orf_id": "te2|orf_1", "domain_types": ""},
        ]
    )
    result = apply_domain_filter(
        df,
        required_domains=["RT"],
        min_domain_types_per_orf=1,
        allow_split_across_orfs=False,
        te_coverage_scope="all",
    )
    row = result.set_index("orf_id")
    assert row.loc["te1|orf_1", "passes_domain_filter"] == True  # noqa: E712
    assert row.loc["te1|orf_2", "passes_domain_filter"] == False  # noqa: E712
    assert row.loc["te2|orf_1", "passes_domain_filter"] == False  # noqa: E712


def test_per_orf_filter_min_two_required():
    df = _make_df(
        [
            {"te_id": "te1", "orf_id": "te1|orf_1", "domain_types": "RT|EN"},
            {"te_id": "te1", "orf_id": "te1|orf_2", "domain_types": "RT"},
        ]
    )
    result = apply_domain_filter(
        df,
        required_domains=["RT", "EN"],
        min_domain_types_per_orf=2,
        allow_split_across_orfs=False,
        te_coverage_scope="all",
    )
    row = result.set_index("orf_id")
    assert row.loc["te1|orf_1", "passes_domain_filter"] == True  # noqa: E712
    assert row.loc["te1|orf_2", "passes_domain_filter"] == False  # noqa: E712


def test_no_split_te_domain_complete_requires_at_least_one_passing_orf():
    df = _make_df(
        [
            {"te_id": "te1", "orf_id": "te1|orf_1", "domain_types": "RT"},
            {"te_id": "te1", "orf_id": "te1|orf_2", "domain_types": ""},
            {"te_id": "te2", "orf_id": "te2|orf_1", "domain_types": ""},
        ]
    )
    result = apply_domain_filter(
        df,
        required_domains=["RT"],
        min_domain_types_per_orf=1,
        allow_split_across_orfs=False,
        te_coverage_scope="all",
    )
    row = result.set_index("orf_id")
    # te1 has one passing ORF → both rows for te1 are complete
    assert row.loc["te1|orf_1", "te_domain_complete"] == True  # noqa: E712
    assert row.loc["te1|orf_2", "te_domain_complete"] == True  # noqa: E712
    # te2 has no passing ORF
    assert row.loc["te2|orf_1", "te_domain_complete"] == False  # noqa: E712


def test_split_across_orfs_all_scope():
    """With split enabled, TE passes when union of all ORFs' domain_types covers required."""
    df = _make_df(
        [
            {
                "te_id": "te1",
                "orf_id": "te1|orf_1",
                "domain_types": "RT",
            },  # passes per-ORF (min=1)
            {
                "te_id": "te1",
                "orf_id": "te1|orf_2",
                "domain_types": "EN",
            },  # passes per-ORF (min=1)
            {"te_id": "te2", "orf_id": "te2|orf_1", "domain_types": "RT"},  # missing EN
        ]
    )
    result = apply_domain_filter(
        df,
        required_domains=["RT", "EN"],
        min_domain_types_per_orf=1,
        allow_split_across_orfs=True,
        te_coverage_scope="all",
    )
    row = result.set_index("orf_id")
    assert row.loc["te1|orf_1", "te_domain_complete"] == True  # noqa: E712
    assert row.loc["te1|orf_2", "te_domain_complete"] == True  # noqa: E712
    assert row.loc["te2|orf_1", "te_domain_complete"] == False  # noqa: E712


def test_split_across_orfs_passing_scope():
    """With te_coverage_scope='passing', only ORFs that pass per-ORF filter contribute to TE coverage."""
    df = _make_df(
        [
            # te1: orf_1 passes (has RT), orf_2 doesn't pass (no required domain) but has EN
            {"te_id": "te1", "orf_id": "te1|orf_1", "domain_types": "RT"},
            {"te_id": "te1", "orf_id": "te1|orf_2", "domain_types": "EN"},
            # te2: both ORFs pass — one has RT, one has EN
            {"te_id": "te2", "orf_id": "te2|orf_1", "domain_types": "RT"},
            {"te_id": "te2", "orf_id": "te2|orf_2", "domain_types": "EN"},
        ]
    )
    # min_domain_types_per_orf=2 so orf_2 of te1 (only EN, not RT) fails per-ORF
    result = apply_domain_filter(
        df,
        required_domains=["RT", "EN"],
        min_domain_types_per_orf=1,
        allow_split_across_orfs=True,
        te_coverage_scope="passing",
    )
    row = result.set_index("orf_id")
    # te1: both orf_1 (RT) and orf_2 (EN) each pass per-ORF (min=1 of {RT,EN}) →
    # union across passing ORFs = {RT, EN} → te_domain_complete=True
    assert row.loc["te1|orf_1", "te_domain_complete"] == True  # noqa: E712
    assert row.loc["te1|orf_2", "te_domain_complete"] == True  # noqa: E712
    # te2: same — both ORFs pass individually and union covers {RT, EN}
    assert row.loc["te2|orf_1", "te_domain_complete"] == True  # noqa: E712
    assert row.loc["te2|orf_2", "te_domain_complete"] == True  # noqa: E712


def test_split_passing_scope_excludes_non_passing_orfs():
    """With te_coverage_scope='passing', an ORF that fails per-ORF does not contribute."""
    df = _make_df(
        [
            # te1: orf_1 has RT (passes), orf_2 has EN but no RT (fails when min=2 of {RT,EN})
            {"te_id": "te1", "orf_id": "te1|orf_1", "domain_types": "RT"},
            {"te_id": "te1", "orf_id": "te1|orf_2", "domain_types": "EN"},
        ]
    )
    result = apply_domain_filter(
        df,
        required_domains=["RT", "EN"],
        min_domain_types_per_orf=2,
        allow_split_across_orfs=True,
        te_coverage_scope="passing",
    )
    row = result.set_index("orf_id")
    # Neither ORF passes per-ORF (each has only 1 of the 2 required types)
    assert row.loc["te1|orf_1", "passes_domain_filter"] == False  # noqa: E712
    assert row.loc["te1|orf_2", "passes_domain_filter"] == False  # noqa: E712
    # No passing ORFs → union is empty → te_domain_complete=False
    assert row.loc["te1|orf_1", "te_domain_complete"] == False  # noqa: E712
    assert row.loc["te1|orf_2", "te_domain_complete"] == False  # noqa: E712
