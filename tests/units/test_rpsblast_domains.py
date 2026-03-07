import os
import sys
from unittest.mock import patch

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
scripts_dir = os.path.join(project_root, "workflow/scripts")
for d in [project_root, scripts_dir]:
    if d not in sys.path:
        sys.path.insert(0, d)

from workflow.scripts.rpsblast_domains import write_empty_outputs, main, RPS_COLS


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _write_aa_fasta(path, ids):
    records = [SeqRecord(Seq("MAAAA"), id=i, description="") for i in ids]
    SeqIO.write(records, str(path), "fasta")


# ---------------------------------------------------------------------------
# write_empty_outputs
# ---------------------------------------------------------------------------

def test_write_empty_outputs_hits_has_header(tmp_path):
    aa = tmp_path / "aa.fa"
    _write_aa_fasta(aa, ["te1|orf_1", "te1|orf_2"])
    hits = str(tmp_path / "hits.tsv")
    summary = str(tmp_path / "summary.tsv")

    write_empty_outputs(str(aa), hits, summary)

    df = pd.read_csv(hits, sep="\t")
    assert list(df.columns) == RPS_COLS + ["qcov"]
    assert df.empty


def test_write_empty_outputs_summary_has_all_orfs(tmp_path):
    aa = tmp_path / "aa.fa"
    _write_aa_fasta(aa, ["te1|orf_1", "te1|orf_2"])
    hits = str(tmp_path / "hits.tsv")
    summary = str(tmp_path / "summary.tsv")

    write_empty_outputs(str(aa), hits, summary)

    df = pd.read_csv(summary, sep="\t")
    assert set(df["orf_id"]) == {"te1|orf_1", "te1|orf_2"}
    assert (df["domain_support"] == False).all()  # noqa: E712
    assert (df["hit_count"] == 0).all()


# ---------------------------------------------------------------------------
# main — no db configured
# ---------------------------------------------------------------------------

def test_main_no_db_writes_empty_outputs(tmp_path):
    aa = tmp_path / "aa.fa"
    _write_aa_fasta(aa, ["te1|orf_1"])
    hits = str(tmp_path / "hits.tsv")
    summary = str(tmp_path / "summary.tsv")

    main(str(aa), hits, summary, db="", evalue=1e-5, max_target_seqs=10,
         min_query_coverage=0.35, min_bitscore=50)

    df = pd.read_csv(summary, sep="\t")
    assert df.iloc[0]["domain_support"] == False  # noqa: E712


# ---------------------------------------------------------------------------
# main — db path not found
# ---------------------------------------------------------------------------

def test_main_db_not_found_writes_empty_outputs(tmp_path):
    aa = tmp_path / "aa.fa"
    _write_aa_fasta(aa, ["te1|orf_1"])
    hits = str(tmp_path / "hits.tsv")
    summary = str(tmp_path / "summary.tsv")
    # use a relative name so Path().glob() doesn't choke on absolute patterns (Python 3.13+)
    main(str(aa), hits, summary, db="nonexistent_db_xyz", evalue=1e-5, max_target_seqs=10,
         min_query_coverage=0.35, min_bitscore=50)

    df = pd.read_csv(summary, sep="\t")
    assert len(df) == 1
    assert df.iloc[0]["domain_support"] == False  # noqa: E712


# ---------------------------------------------------------------------------
# main — db found, subprocess mocked
# ---------------------------------------------------------------------------

RAW_HITS = "\t".join([
    "te1|orf_1", "CDD:12345", "95.0", "150", "5", "0",
    "1", "150", "1", "150", "1e-50", "300.0", "Reverse transcriptase",
]) + "\n"


def _fake_rpsblast_empty(raw_out):
    """Mock: write an empty hits file."""
    open(raw_out, "w").close()


def _fake_rpsblast_with_hits(raw_out):
    """Mock: write one hit."""
    with open(raw_out, "w") as fh:
        fh.write(RAW_HITS)


def test_main_empty_blast_output(tmp_path):
    aa = tmp_path / "aa.fa"
    _write_aa_fasta(aa, ["te1|orf_1"])
    db_file = tmp_path / "mydb"
    db_file.write_text("")  # make it exist
    hits = str(tmp_path / "hits.tsv")
    summary = str(tmp_path / "summary.tsv")

    def fake_run(cmd, check):
        raw_out = cmd[cmd.index("-out") + 1]
        _fake_rpsblast_empty(raw_out)

    with patch("workflow.scripts.rpsblast_domains.subprocess.run", side_effect=fake_run):
        main(str(aa), hits, summary, db=str(db_file), evalue=1e-5, max_target_seqs=10,
             min_query_coverage=0.35, min_bitscore=50)

    df = pd.read_csv(summary, sep="\t")
    assert df.iloc[0]["domain_support"] == False  # noqa: E712


def test_main_with_hits_sets_domain_support(tmp_path):
    aa = tmp_path / "aa.fa"
    _write_aa_fasta(aa, ["te1|orf_1"])
    db_file = tmp_path / "mydb"
    db_file.write_text("")
    hits = str(tmp_path / "hits.tsv")
    summary = str(tmp_path / "summary.tsv")

    def fake_run(cmd, check):
        raw_out = cmd[cmd.index("-out") + 1]
        _fake_rpsblast_with_hits(raw_out)

    with patch("workflow.scripts.rpsblast_domains.subprocess.run", side_effect=fake_run):
        main(str(aa), hits, summary, db=str(db_file), evalue=1e-5, max_target_seqs=10,
             min_query_coverage=0.35, min_bitscore=50)

    df = pd.read_csv(summary, sep="\t")
    row = df[df["orf_id"] == "te1|orf_1"].iloc[0]
    assert row["domain_support"] == True  # noqa: E712
    assert row["best_domain"] == "Reverse transcriptase"
    assert int(row["hit_count"]) == 1


def test_main_hit_filtered_by_bitscore(tmp_path):
    aa = tmp_path / "aa.fa"
    _write_aa_fasta(aa, ["te1|orf_1"])
    db_file = tmp_path / "mydb"
    db_file.write_text("")
    hits = str(tmp_path / "hits.tsv")
    summary = str(tmp_path / "summary.tsv")

    low_bitscore_hit = "\t".join([
        "te1|orf_1", "CDD:12345", "95.0", "150", "5", "0",
        "1", "150", "1", "150", "1e-50", "10.0", "Reverse transcriptase",
    ]) + "\n"

    def fake_run(cmd, check):
        raw_out = cmd[cmd.index("-out") + 1]
        with open(raw_out, "w") as fh:
            fh.write(low_bitscore_hit)

    with patch("workflow.scripts.rpsblast_domains.subprocess.run", side_effect=fake_run):
        main(str(aa), hits, summary, db=str(db_file), evalue=1e-5, max_target_seqs=10,
             min_query_coverage=0.35, min_bitscore=50)

    df = pd.read_csv(summary, sep="\t")
    assert df.iloc[0]["domain_support"] == False  # noqa: E712
