import os
import sys
from io import StringIO
from Bio import SeqIO
import pytest

# Add the project root to sys.path to allow importing from workflow.scripts
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from workflow.scripts.collect_te import collect_te, main


@pytest.fixture
def mock_fasta_content():
    return (
        ">TE1_sequence1\n"
        "ATGC\n"
        ">TE1_sequence2\n"
        "GCAT\n"
        ">TE2_sequence1\n"
        "AAAA\n"
        ">OTHER_sequence1\n"
        "TTTT\n"
    )


@pytest.fixture
def mock_fasta_file(mock_fasta_content):
    return StringIO(mock_fasta_content)


def test_collect_te_found(mock_fasta_file):
    results = collect_te(mock_fasta_file, "TE1")
    assert len(results) == 2
    assert results[0].id == "TE1_sequence1"
    assert results[1].id == "TE1_sequence2"
    assert str(results[0].seq) == "ATGC"


def test_collect_te_not_found(mock_fasta_file):
    results = collect_te(mock_fasta_file, "TE3")
    assert len(results) == 0


def test_collect_te_partial_match(mock_fasta_file):
    results = collect_te(mock_fasta_file, "TE")
    assert len(results) == 3
    ids = [r.id for r in results]
    assert "TE1_sequence1" in ids
    assert "TE1_sequence2" in ids
    assert "TE2_sequence1" in ids


def test_collect_te_empty_fasta():
    empty_fasta = StringIO("")
    results = collect_te(empty_fasta, "TE1")
    assert len(results) == 0


def test_main_logic(tmp_path, mock_fasta_content):
    # Setup temporary files
    input_file = tmp_path / "input.fasta"
    output_file = tmp_path / "output.fasta"

    input_file.write_text(mock_fasta_content)

    # Run main logic
    main(str(input_file), str(output_file), "TE1")

    # Verify output
    assert output_file.exists()
    records = list(SeqIO.parse(str(output_file), "fasta"))
    assert len(records) == 2
    assert records[0].id == "TE1_sequence1"
    assert records[1].id == "TE1_sequence2"
