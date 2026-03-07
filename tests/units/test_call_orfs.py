import os
import sys
from io import StringIO
from unittest.mock import call, patch

import pytest
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
scripts_dir = os.path.join(project_root, "workflow/scripts")
for d in [project_root, scripts_dir]:
    if d not in sys.path:
        sys.path.insert(0, d)

from workflow.scripts.call_orfs import parse_getorf_header, run_getorf, main


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def make_record(description):
    """Build a minimal SeqRecord whose .description matches getorf output."""
    first_word = description.split()[0]
    return SeqRecord(Seq("ATGCCC"), id=first_word, description=description)


# ---------------------------------------------------------------------------
# parse_getorf_header
# ---------------------------------------------------------------------------

def test_parse_forward_strand():
    rec = make_record("te1_1 [10 - 100]")
    result = parse_getorf_header(rec)
    assert result["te_id"] == "te1"
    assert result["start"] == 10
    assert result["end"] == 100
    assert result["strand"] == "+"
    assert result["frame"] is None


def test_parse_reverse_strand():
    rec = make_record("te1_2 [200 - 100] (REVERSE SENSE)")
    result = parse_getorf_header(rec)
    assert result["te_id"] == "te1"
    assert result["start"] == 200
    assert result["end"] == 100
    assert result["strand"] == "-"


def test_parse_seqid_with_underscores():
    """Seqid containing underscores: only the trailing _N suffix must be stripped."""
    rec = make_record("my_te_seq_3 [50 - 150]")
    result = parse_getorf_header(rec)
    assert result["te_id"] == "my_te_seq"


def test_parse_missing_coordinates():
    rec = make_record("te1_1 no coords here")
    result = parse_getorf_header(rec)
    assert result["start"] is None
    assert result["end"] is None


def test_parse_reverse_sense_case_insensitive():
    rec = make_record("te1_1 [10 - 90] (reverse sense)")
    assert parse_getorf_header(rec)["strand"] == "-"


# ---------------------------------------------------------------------------
# run_getorf
# ---------------------------------------------------------------------------

def test_run_getorf_canonical(tmp_path):
    out = str(tmp_path / "out.fa")
    with patch("workflow.scripts.call_orfs.subprocess.run") as mock_run:
        run_getorf("input.fa", out, find=3, min_orf_nt=300)
    mock_run.assert_called_once_with(
        ["getorf", "-sequence", "input.fa", "-outseq", out, "-minsize", "300", "-find", "3"],
        check=True,
    )


def test_run_getorf_permissive(tmp_path):
    out = str(tmp_path / "out.fa")
    with patch("workflow.scripts.call_orfs.subprocess.run") as mock_run:
        run_getorf("input.fa", out, find=0, min_orf_nt=150)
    args = mock_run.call_args[0][0]
    assert "-find" in args
    assert args[args.index("-find") + 1] == "0"


def test_run_getorf_propagates_failure(tmp_path):
    import subprocess
    out = str(tmp_path / "out.fa")
    with patch("workflow.scripts.call_orfs.subprocess.run", side_effect=subprocess.CalledProcessError(1, "getorf")):
        with pytest.raises(subprocess.CalledProcessError):
            run_getorf("input.fa", out, find=1, min_orf_nt=300)


# ---------------------------------------------------------------------------
# main — helpers to fake getorf output
# ---------------------------------------------------------------------------

NT_FASTA = """\
>te1_1 [10 - 100]
ATGCCCAAA
>te1_2 [200 - 110] (REVERSE SENSE)
ATGAAACCC
"""

AA_FASTA = """\
>te1_1 [10 - 100]
MPA
>te1_2 [200 - 110] (REVERSE SENSE)
MKP
"""


def fake_run_getorf(input_fasta, out_path, find, min_orf_nt):
    """Write pre-baked FASTA so main() has something to parse."""
    content = NT_FASTA if find in (2, 3) else AA_FASTA
    with open(out_path, "w") as fh:
        fh.write(content)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def test_main_canonical_mode(tmp_path):
    table = str(tmp_path / "orfs.tsv")
    nt_fa = str(tmp_path / "nt.fa")
    aa_fa = str(tmp_path / "aa.fa")

    with patch("workflow.scripts.call_orfs.run_getorf", side_effect=fake_run_getorf) as mock_getorf:
        main("input.fa", table, nt_fa, aa_fa, min_orf_aa=3, require_start_codon=True)

    # canonical mode → find=3 (nt) and find=1 (aa)
    calls = mock_getorf.call_args_list
    assert calls[0].kwargs["find"] == 3
    assert calls[1].kwargs["find"] == 1

    import pandas as pd
    df = pd.read_csv(table, sep="\t")
    assert len(df) == 2
    assert df["has_start"].all()


def test_main_permissive_mode(tmp_path):
    table = str(tmp_path / "orfs.tsv")
    nt_fa = str(tmp_path / "nt.fa")
    aa_fa = str(tmp_path / "aa.fa")

    with patch("workflow.scripts.call_orfs.run_getorf", side_effect=fake_run_getorf) as mock_getorf:
        main("input.fa", table, nt_fa, aa_fa, min_orf_aa=3, require_start_codon=False)

    calls = mock_getorf.call_args_list
    assert calls[0].kwargs["find"] == 2
    assert calls[1].kwargs["find"] == 0

    import pandas as pd
    df = pd.read_csv(table, sep="\t")
    assert not df["has_start"].any()


def test_main_output_columns(tmp_path):
    table = str(tmp_path / "orfs.tsv")
    with patch("workflow.scripts.call_orfs.run_getorf", side_effect=fake_run_getorf):
        main("input.fa", table, str(tmp_path / "nt.fa"), str(tmp_path / "aa.fa"), min_orf_aa=3)

    import pandas as pd
    df = pd.read_csv(table, sep="\t")
    for col in ("te_id", "orf_id", "start", "end", "strand", "frame", "nt_len", "aa_len", "has_start", "has_stop"):
        assert col in df.columns


def test_main_strand_assignment(tmp_path):
    table = str(tmp_path / "orfs.tsv")
    with patch("workflow.scripts.call_orfs.run_getorf", side_effect=fake_run_getorf):
        main("input.fa", table, str(tmp_path / "nt.fa"), str(tmp_path / "aa.fa"), min_orf_aa=3)

    import pandas as pd
    df = pd.read_csv(table, sep="\t")
    assert df.iloc[0]["strand"] == "+"
    assert df.iloc[1]["strand"] == "-"


def test_main_empty_input(tmp_path):
    """When getorf produces no ORFs, outputs should be empty but valid."""
    table = str(tmp_path / "orfs.tsv")
    nt_fa = str(tmp_path / "nt.fa")
    aa_fa = str(tmp_path / "aa.fa")

    def empty_getorf(input_fasta, out_path, find, min_orf_nt):
        open(out_path, "w").close()

    with patch("workflow.scripts.call_orfs.run_getorf", side_effect=empty_getorf):
        main("input.fa", table, nt_fa, aa_fa, min_orf_aa=3)

    import pandas as pd
    df = pd.read_csv(table, sep="\t")
    assert df.empty
    for col in ("te_id", "orf_id", "start", "end", "strand", "frame", "nt_len", "aa_len", "has_start", "has_stop"):
        assert col in df.columns


def test_main_count_mismatch_logs_warning(tmp_path, caplog):
    """A mismatch between nt and aa record counts should trigger a warning."""
    import logging

    def mismatched_getorf(input_fasta, out_path, find, min_orf_nt):
        if find in (2, 3):
            with open(out_path, "w") as fh:
                fh.write(NT_FASTA)
        else:
            with open(out_path, "w") as fh:
                fh.write(AA_FASTA.split("\n>")[0] + "\n")  # only first record

    with patch("workflow.scripts.call_orfs.run_getorf", side_effect=mismatched_getorf):
        with caplog.at_level(logging.WARNING, logger="workflow.scripts.call_orfs"):
            main("input.fa", str(tmp_path / "t.tsv"), str(tmp_path / "nt.fa"), str(tmp_path / "aa.fa"), min_orf_aa=1)

    assert any("mismatch" in r.message.lower() for r in caplog.records)
