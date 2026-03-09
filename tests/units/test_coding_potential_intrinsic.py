import os
import sys

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

from workflow.scripts.coding_potential_intrinsic import gc_content, codon_diversity, intrinsic_label, main


# ---------------------------------------------------------------------------
# gc_content
# ---------------------------------------------------------------------------

def test_gc_content_all_gc():
    assert gc_content("GCGCGC") == 1.0


def test_gc_content_all_at():
    assert gc_content("ATATAT") == 0.0


def test_gc_content_mixed():
    assert gc_content("ATGC") == pytest.approx(0.5)


def test_gc_content_empty():
    assert gc_content("") == 0.0


def test_gc_content_case_insensitive():
    assert gc_content("atgc") == pytest.approx(0.5)


# ---------------------------------------------------------------------------
# codon_diversity
# ---------------------------------------------------------------------------

def test_codon_diversity_all_same():
    # ATG ATG ATG → 1 unique / 3 total
    assert codon_diversity("ATGATGATG") == pytest.approx(1 / 3)


def test_codon_diversity_all_different():
    # ATG CCC AAA → 3 unique / 3 total
    assert codon_diversity("ATGCCCAAA") == pytest.approx(1.0)


def test_codon_diversity_empty():
    assert codon_diversity("") == 0.0


def test_codon_diversity_too_short():
    assert codon_diversity("AT") == 0.0


def test_codon_diversity_partial_trailing_ignored():
    # 3 complete codons + 1 leftover nucleotide
    assert codon_diversity("ATGCCCAAAT") == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# intrinsic_label
# ---------------------------------------------------------------------------

def test_label_high():
    assert intrinsic_label(0.8, high_threshold=0.7, medium_threshold=0.45) == "high"


def test_label_medium():
    assert intrinsic_label(0.6, high_threshold=0.7, medium_threshold=0.45) == "medium"


def test_label_low():
    assert intrinsic_label(0.3, high_threshold=0.7, medium_threshold=0.45) == "low"


def test_label_boundary_high():
    assert intrinsic_label(0.7, high_threshold=0.7, medium_threshold=0.45) == "high"


def test_label_boundary_medium():
    assert intrinsic_label(0.45, high_threshold=0.7, medium_threshold=0.45) == "medium"


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def _make_inputs(tmp_path, rows, fasta_seqs):
    table = tmp_path / "orfs.tsv"
    fasta = tmp_path / "orfs.fa"
    pd.DataFrame(rows).to_csv(str(table), sep="\t", index=False)
    records = [SeqRecord(Seq(seq), id=oid, description="") for oid, seq in fasta_seqs.items()]
    SeqIO.write(records, str(fasta), "fasta")
    return str(table), str(fasta)


def test_main_output_columns(tmp_path):
    rows = [{"orf_id": "te1|orf_1", "aa_len": 150}]
    seqs = {"te1|orf_1": "ATGCCCAAAGGG"}
    table, fasta = _make_inputs(tmp_path, rows, seqs)
    out = str(tmp_path / "scores.tsv")

    main(table, fasta, out, high_threshold=0.7, medium_threshold=0.45)

    df = pd.read_csv(out, sep="\t")
    for col in ("orf_id", "gc_content", "codon_diversity", "length_score", "intrinsic_score", "intrinsic_label"):
        assert col in df.columns


def test_main_score_values_in_range(tmp_path):
    rows = [{"orf_id": "te1|orf_1", "aa_len": 300}]
    seqs = {"te1|orf_1": "ATGCCCGGG" * 10}
    table, fasta = _make_inputs(tmp_path, rows, seqs)
    out = str(tmp_path / "scores.tsv")

    main(table, fasta, out, high_threshold=0.7, medium_threshold=0.45)

    df = pd.read_csv(out, sep="\t")
    assert 0.0 <= df.iloc[0]["intrinsic_score"] <= 1.0
    assert 0.0 <= df.iloc[0]["gc_content"] <= 1.0
    assert df.iloc[0]["intrinsic_label"] in ("high", "medium", "low")


def test_main_length_score_capped_at_one(tmp_path):
    # aa_len=600 → length_score = min(600/300, 1.0) = 1.0
    rows = [{"orf_id": "te1|orf_1", "aa_len": 600}]
    seqs = {"te1|orf_1": "ATG" * 200}
    table, fasta = _make_inputs(tmp_path, rows, seqs)
    out = str(tmp_path / "scores.tsv")

    main(table, fasta, out, high_threshold=0.7, medium_threshold=0.45)

    df = pd.read_csv(out, sep="\t")
    assert df.iloc[0]["length_score"] == pytest.approx(1.0)


def test_main_missing_fasta_id_raises(tmp_path):
    rows = [{"orf_id": "te1|orf_1", "aa_len": 150}]
    seqs = {}  # orf_id absent from FASTA
    table, fasta = _make_inputs(tmp_path, rows, seqs)

    with pytest.raises(ValueError, match="missing from FASTA"):
        main(table, fasta, str(tmp_path / "scores.tsv"), high_threshold=0.7, medium_threshold=0.45)


def test_main_invalid_threshold_out_of_range(tmp_path):
    rows = [{"orf_id": "te1|orf_1", "aa_len": 150}]
    seqs = {"te1|orf_1": "ATGCCC"}
    table, fasta = _make_inputs(tmp_path, rows, seqs)
    with pytest.raises(ValueError, match="Thresholds must be in"):
        main(table, fasta, str(tmp_path / "scores.tsv"), high_threshold=1.5, medium_threshold=0.45)


def test_main_invalid_threshold_medium_gt_high(tmp_path):
    rows = [{"orf_id": "te1|orf_1", "aa_len": 150}]
    seqs = {"te1|orf_1": "ATGCCC"}
    table, fasta = _make_inputs(tmp_path, rows, seqs)
    with pytest.raises(ValueError, match="medium_threshold"):
        main(table, fasta, str(tmp_path / "scores.tsv"), high_threshold=0.4, medium_threshold=0.7)


def test_main_empty_table(tmp_path):
    # write a table with column headers but no rows
    table = str(tmp_path / "orfs.tsv")
    fasta = str(tmp_path / "orfs.fa")
    pd.DataFrame(columns=["orf_id", "aa_len"]).to_csv(table, sep="\t", index=False)
    open(fasta, "w").close()
    out = str(tmp_path / "scores.tsv")

    main(table, fasta, out, high_threshold=0.7, medium_threshold=0.45)

    df = pd.read_csv(out, sep="\t")
    assert df.empty
    for col in ("orf_id", "gc_content", "codon_diversity", "length_score", "intrinsic_score", "intrinsic_label"):
        assert col in df.columns
