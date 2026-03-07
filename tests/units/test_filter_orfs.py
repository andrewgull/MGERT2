import os
import sys
from io import StringIO

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

from workflow.scripts.filter_orfs import load_fasta_dict, write_subset_fasta, main


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def fasta_file(tmp_path):
    path = tmp_path / "seqs.fa"
    records = [
        SeqRecord(Seq("ATGCCCAAA"), id="te1|orf_1", description=""),
        SeqRecord(Seq("ATGAAACCC"), id="te1|orf_2", description=""),
        SeqRecord(Seq("GCTAGCTA"), id="te2|orf_1", description=""),
    ]
    SeqIO.write(records, str(path), "fasta")
    return path


@pytest.fixture
def orf_table(tmp_path):
    df = pd.DataFrame(
        {
            "te_id": ["te1", "te1", "te1", "te2"],
            "orf_id": ["te1|orf_1", "te1|orf_2", "te1|orf_3", "te2|orf_1"],
            "aa_len": [200, 150, 80, 120],
            "has_stop": [True, True, True, True],
        }
    )
    path = tmp_path / "orfs.tsv"
    df.to_csv(str(path), sep="\t", index=False)
    return path


# ---------------------------------------------------------------------------
# load_fasta_dict
# ---------------------------------------------------------------------------

def test_load_fasta_dict(fasta_file):
    d = load_fasta_dict(str(fasta_file))
    assert set(d.keys()) == {"te1|orf_1", "te1|orf_2", "te2|orf_1"}
    assert str(d["te1|orf_1"].seq) == "ATGCCCAAA"


# ---------------------------------------------------------------------------
# write_subset_fasta
# ---------------------------------------------------------------------------

def test_write_subset_fasta(tmp_path, fasta_file):
    records = load_fasta_dict(str(fasta_file))
    out = str(tmp_path / "out.fa")
    write_subset_fasta(records, {"te1|orf_1"}, out)
    written = list(SeqIO.parse(out, "fasta"))
    assert len(written) == 1
    assert written[0].id == "te1|orf_1"


def test_write_subset_fasta_missing_id_skipped(tmp_path, fasta_file):
    records = load_fasta_dict(str(fasta_file))
    out = str(tmp_path / "out.fa")
    write_subset_fasta(records, {"te1|orf_1", "nonexistent"}, out)
    written = list(SeqIO.parse(out, "fasta"))
    assert len(written) == 1


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def _write_fasta(path, ids):
    records = [SeqRecord(Seq("ATGCCC"), id=i, description="") for i in ids]
    SeqIO.write(records, str(path), "fasta")


def _make_table(tmp_path, rows):
    path = tmp_path / "in.tsv"
    pd.DataFrame(rows).to_csv(str(path), sep="\t", index=False)
    return path


def test_main_filters_by_aa_len(tmp_path):
    rows = [
        {"te_id": "te1", "orf_id": "te1|orf_1", "aa_len": 200, "has_stop": True},
        {"te_id": "te1", "orf_id": "te1|orf_2", "aa_len": 50, "has_stop": True},
    ]
    table = _make_table(tmp_path, rows)
    nt = tmp_path / "nt.fa"
    aa = tmp_path / "aa.fa"
    _write_fasta(nt, ["te1|orf_1", "te1|orf_2"])
    _write_fasta(aa, ["te1|orf_1", "te1|orf_2"])
    out_table = str(tmp_path / "out.tsv")

    main(str(table), str(nt), str(aa), out_table, str(tmp_path / "nt_out.fa"), str(tmp_path / "aa_out.fa"),
         min_orf_aa=100, max_orfs_per_te=5)

    df = pd.read_csv(out_table, sep="\t")
    assert list(df["orf_id"]) == ["te1|orf_1"]


def test_main_require_stop_codon_filters(tmp_path):
    rows = [
        {"te_id": "te1", "orf_id": "te1|orf_1", "aa_len": 200, "has_stop": True},
        {"te_id": "te1", "orf_id": "te1|orf_2", "aa_len": 200, "has_stop": False},
    ]
    table = _make_table(tmp_path, rows)
    nt = tmp_path / "nt.fa"
    aa = tmp_path / "aa.fa"
    _write_fasta(nt, ["te1|orf_1", "te1|orf_2"])
    _write_fasta(aa, ["te1|orf_1", "te1|orf_2"])
    out_table = str(tmp_path / "out.tsv")

    main(str(table), str(nt), str(aa), out_table, str(tmp_path / "nt_out.fa"), str(tmp_path / "aa_out.fa"),
         min_orf_aa=100, max_orfs_per_te=5, require_stop_codon=True)

    df = pd.read_csv(out_table, sep="\t")
    assert list(df["orf_id"]) == ["te1|orf_1"]


def test_main_no_require_stop_codon_keeps_all(tmp_path):
    rows = [
        {"te_id": "te1", "orf_id": "te1|orf_1", "aa_len": 200, "has_stop": False},
    ]
    table = _make_table(tmp_path, rows)
    nt = tmp_path / "nt.fa"
    aa = tmp_path / "aa.fa"
    _write_fasta(nt, ["te1|orf_1"])
    _write_fasta(aa, ["te1|orf_1"])
    out_table = str(tmp_path / "out.tsv")

    main(str(table), str(nt), str(aa), out_table, str(tmp_path / "nt_out.fa"), str(tmp_path / "aa_out.fa"),
         min_orf_aa=100, max_orfs_per_te=5, require_stop_codon=False)

    df = pd.read_csv(out_table, sep="\t")
    assert len(df) == 1


def test_main_max_orfs_per_te(tmp_path):
    rows = [
        {"te_id": "te1", "orf_id": f"te1|orf_{i}", "aa_len": 300 - i * 10, "has_stop": True}
        for i in range(5)
    ]
    table = _make_table(tmp_path, rows)
    ids = [r["orf_id"] for r in rows]
    nt = tmp_path / "nt.fa"
    aa = tmp_path / "aa.fa"
    _write_fasta(nt, ids)
    _write_fasta(aa, ids)
    out_table = str(tmp_path / "out.tsv")

    main(str(table), str(nt), str(aa), out_table, str(tmp_path / "nt_out.fa"), str(tmp_path / "aa_out.fa"),
         min_orf_aa=100, max_orfs_per_te=2)

    df = pd.read_csv(out_table, sep="\t")
    assert len(df) == 2
    # top 2 by aa_len
    assert df.iloc[0]["aa_len"] >= df.iloc[1]["aa_len"]


def test_main_empty_input(tmp_path):
    # write a table with column headers but no rows so pandas can read it
    table = tmp_path / "in.tsv"
    pd.DataFrame(columns=["te_id", "orf_id", "aa_len", "has_stop"]).to_csv(str(table), sep="\t", index=False)
    nt = tmp_path / "nt.fa"
    aa = tmp_path / "aa.fa"
    nt.write_text("")
    aa.write_text("")
    out_table = str(tmp_path / "out.tsv")
    out_nt = str(tmp_path / "nt_out.fa")
    out_aa = str(tmp_path / "aa_out.fa")

    main(str(table), str(nt), str(aa), out_table, out_nt, out_aa, min_orf_aa=100, max_orfs_per_te=3)

    df = pd.read_csv(out_table, sep="\t")
    assert df.empty


def test_main_creates_output_dirs(tmp_path):
    """Parent directories for all outputs are created even when they don't exist."""
    rows = [{"te_id": "te1", "orf_id": "te1|orf_1", "aa_len": 200, "has_stop": True}]
    table = _make_table(tmp_path, rows)
    nested = tmp_path / "a" / "b" / "c"
    nt = nested / "nt.fa"
    aa = nested / "aa.fa"
    _write_fasta(tmp_path / "nt_in.fa", ["te1|orf_1"])
    _write_fasta(tmp_path / "aa_in.fa", ["te1|orf_1"])
    out_table = str(nested / "out.tsv")

    main(str(table), str(tmp_path / "nt_in.fa"), str(tmp_path / "aa_in.fa"),
         out_table, str(nt), str(aa), min_orf_aa=100, max_orfs_per_te=5)

    assert (nested / "out.tsv").exists()


def test_main_creates_output_dirs_on_empty_input(tmp_path):
    """mkdir is called even in the early-return (empty input) branch."""
    table = tmp_path / "in.tsv"
    pd.DataFrame(columns=["te_id", "orf_id", "aa_len", "has_stop"]).to_csv(str(table), sep="\t", index=False)
    nested = tmp_path / "x" / "y"

    main(str(table), str(tmp_path / "nt.fa"), str(tmp_path / "aa.fa"),
         str(nested / "out.tsv"), str(nested / "nt_out.fa"), str(nested / "aa_out.fa"),
         min_orf_aa=100, max_orfs_per_te=3)

    assert (nested / "out.tsv").exists()


def test_main_fasta_outputs_match_filtered(tmp_path):
    rows = [
        {"te_id": "te1", "orf_id": "te1|orf_1", "aa_len": 200, "has_stop": True},
        {"te_id": "te1", "orf_id": "te1|orf_2", "aa_len": 50, "has_stop": True},
    ]
    table = _make_table(tmp_path, rows)
    nt = tmp_path / "nt.fa"
    aa = tmp_path / "aa.fa"
    _write_fasta(nt, ["te1|orf_1", "te1|orf_2"])
    _write_fasta(aa, ["te1|orf_1", "te1|orf_2"])
    out_nt = str(tmp_path / "nt_out.fa")
    out_aa = str(tmp_path / "aa_out.fa")

    main(str(table), str(nt), str(aa), str(tmp_path / "out.tsv"), out_nt, out_aa,
         min_orf_aa=100, max_orfs_per_te=5)

    nt_ids = [r.id for r in SeqIO.parse(out_nt, "fasta")]
    aa_ids = [r.id for r in SeqIO.parse(out_aa, "fasta")]
    assert nt_ids == ["te1|orf_1"]
    assert aa_ids == ["te1|orf_1"]
