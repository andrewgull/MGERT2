import pytest
import pandas as pd
import numpy as np
import os
import gzip
from workflow.scripts.repeat_landscape import (
    get_genome_size,
    parse_repeatmasker_out,
    plot_landscape,
)

@pytest.fixture
def mock_fasta(tmp_path):
    """Creates a mock FASTA file."""
    fasta_path = tmp_path / "test.fasta"
    content = ">seq1\nATGC\n>seq2\nNNNN\n"
    fasta_path.write_text(content)
    return str(fasta_path)

@pytest.fixture
def mock_fasta_gz(tmp_path):
    """Creates a mock gzipped FASTA file."""
    fasta_path = tmp_path / "test.fasta.gz"
    content = ">seq1\nATGC\n>seq2\nNNNN\n".encode("utf-8")
    with gzip.open(fasta_path, "wb") as f:
        f.write(content)
    return str(fasta_path)

@pytest.fixture
def mock_rm_out(tmp_path):
    """Creates a mock RepeatMasker .out file."""
    out_path = tmp_path / "test.out"
    # RepeatMasker .out format is fixed-width. 
    # Header lines (3) then columns.
    header = [
        "   SW   perc perc perc  query      position in query           matching          repeat              position in  repeat\n",
        "score   div. del. ins.  sequence    begin     end    (left)    repeat            class/family         begin  end (left)   ID\n",
        "\n"
    ]
    # Example entries
    # 1. LINE/L1, div 5, length 100
    # 2. SINE/Alu, div 10, length 50
    # cols = ["sw_score", "perc_div", "perc_del", "perc_ins", "query", "q_start", "q_end", "q_left", "strand", "repeat", "class_family", "r_start", "r_end", "r_left", "id"]
    line1 = "  100    5.0  0.0  0.0  chr1           1     100   (1000)  +  L1              LINE/L1               1  100    (0)    1\n"
    line2 = "  200   10.0  0.0  0.0  chr1         200     250    (750)  +  Alu             SINE/Alu              1   50    (0)    2\n"
    
    with open(out_path, "w") as f:
        f.writelines(header)
        f.write(line1)
        f.write(line2)
    
    return str(out_path)

def test_get_genome_size(mock_fasta, mock_fasta_gz):
    # ATGC (4) + NNNN (4) = 8
    assert get_genome_size(mock_fasta) == 8
    assert get_genome_size(mock_fasta_gz) == 8

def test_parse_repeatmasker_out(mock_rm_out):
    genome_size = 1000
    df = parse_repeatmasker_out(mock_rm_out, genome_size)
    
    # Check shape/columns
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    
    # L1: len 99 (100-1), div 5.0 -> bin 5.0
    # Alu: len 50 (250-200), div 10.0 -> bin 10.0
    # genome_perc calculation: (length / 1000) * 100
    
    l1_data = df[df["class_family"] == "LINE/L1"]
    assert len(l1_data) == 1
    # q_end - q_start: 100 - 1 = 99
    assert l1_data.iloc[0]["length"] == 100
    assert l1_data.iloc[0]["div_bin"] == 5.0
    assert l1_data.iloc[0]["genome_perc"] == (100 / 1000) * 100

def test_plot_landscape(tmp_path, mock_rm_out):
    genome_size = 1000
    df = parse_repeatmasker_out(mock_rm_out, genome_size)
    output_img = tmp_path / "landscape.png"
    
    # Just check it runs and saves
    plot_landscape(df, str(output_img))
    assert output_img.exists()

def test_main_integration(tmp_path, mock_fasta, mock_rm_out):
    # Import main from the script
    from workflow.scripts.repeat_landscape import main
    
    output_img = tmp_path / "main_landscape.png"
    # Now main takes a directory (tmp_path contains mock_rm_out)
    main(mock_fasta, str(tmp_path), str(output_img))
    
    assert output_img.exists()

def test_main_no_out_file(tmp_path, mock_fasta):
    from workflow.scripts.repeat_landscape import main
    import pytest
    
    output_img = tmp_path / "fail_landscape.png"
    # empty_dir has no .out files
    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()
    
    with pytest.raises(FileNotFoundError, match="No .out file found"):
        main(mock_fasta, str(empty_dir), str(output_img))
