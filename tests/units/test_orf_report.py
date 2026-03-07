import os
import sys

import pandas as pd
import pytest

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
scripts_dir = os.path.join(project_root, "workflow/scripts")
for d in [project_root, scripts_dir]:
    if d not in sys.path:
        sys.path.insert(0, d)

from workflow.scripts.orf_report import main


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _write_tsv(path, df):
    df.to_csv(str(path), sep="\t", index=False)


@pytest.fixture
def report_inputs(tmp_path):
    raw = pd.DataFrame({
        "te_id": ["te1", "te1", "te2"],
        "orf_id": ["te1|orf_1", "te1|orf_2", "te2|orf_1"],
        "aa_len": [400, 200, 100],
    })
    classified = pd.DataFrame({
        "te_id": ["te1", "te1"],
        "orf_id": ["te1|orf_1", "te1|orf_2"],
        "aa_len": [400, 200],
        "confidence_class": ["high_confidence_coding", "putative_coding"],
    })
    _write_tsv(tmp_path / "raw.tsv", raw)
    _write_tsv(tmp_path / "classified.tsv", classified)
    return tmp_path


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def test_main_summary_tsv_written(report_inputs):
    out_tsv = str(report_inputs / "summary.tsv")
    out_html = str(report_inputs / "summary.html")
    main(str(report_inputs / "raw.tsv"), str(report_inputs / "classified.tsv"), out_tsv, out_html)

    df = pd.read_csv(out_tsv, sep="\t")
    metrics = dict(zip(df["metric"], df["value"]))
    assert int(metrics["total_raw_orfs"]) == 3
    assert int(metrics["total_filtered_orfs"]) == 2
    assert int(metrics["te_with_raw_orfs"]) == 2
    assert int(metrics["te_with_filtered_orfs"]) == 1
    assert int(metrics["longest_orf_aa"]) == 400


def test_main_class_counts_in_summary(report_inputs):
    out_tsv = str(report_inputs / "summary.tsv")
    out_html = str(report_inputs / "summary.html")
    main(str(report_inputs / "raw.tsv"), str(report_inputs / "classified.tsv"), out_tsv, out_html)

    df = pd.read_csv(out_tsv, sep="\t")
    metrics = dict(zip(df["metric"], df["value"]))
    assert int(metrics["high_confidence_coding"]) == 1
    assert int(metrics["putative_coding"]) == 1


def test_main_html_written(report_inputs):
    out_tsv = str(report_inputs / "summary.tsv")
    out_html = str(report_inputs / "summary.html")
    main(str(report_inputs / "raw.tsv"), str(report_inputs / "classified.tsv"), out_tsv, out_html)

    html = open(out_html).read()
    assert "<html>" in html
    assert "ORF Coding Potential Report" in html
    assert "<table" in html


def test_main_empty_classified(tmp_path):
    raw = pd.DataFrame({"te_id": ["te1"], "orf_id": ["te1|orf_1"], "aa_len": [200]})
    classified = pd.DataFrame(columns=["te_id", "orf_id", "aa_len", "confidence_class"])
    _write_tsv(tmp_path / "raw.tsv", raw)
    _write_tsv(tmp_path / "classified.tsv", classified)
    out_tsv = str(tmp_path / "summary.tsv")
    out_html = str(tmp_path / "summary.html")

    main(str(tmp_path / "raw.tsv"), str(tmp_path / "classified.tsv"), out_tsv, out_html)

    df = pd.read_csv(out_tsv, sep="\t")
    metrics = dict(zip(df["metric"], df["value"]))
    assert int(metrics["longest_orf_aa"]) == 0
    assert int(metrics["total_filtered_orfs"]) == 0


def test_main_creates_output_dir(tmp_path, report_inputs):
    out_tsv = str(tmp_path / "nested" / "summary.tsv")
    out_html = str(tmp_path / "nested" / "summary.html")
    main(str(report_inputs / "raw.tsv"), str(report_inputs / "classified.tsv"), out_tsv, out_html)
    assert os.path.exists(out_tsv)
