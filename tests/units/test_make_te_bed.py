import os
import sys
import pandas as pd
import pytest
from io import StringIO

# Add the project root and scripts directory to sys.path to allow imports
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
scripts_dir = os.path.join(project_root, "workflow/scripts")
for d in [project_root, scripts_dir]:
    if d not in sys.path:
        sys.path.insert(0, d)

from workflow.scripts.make_te_bed import parse_repeatmasker_out, filter_by_te_name, convert_to_bed, write_bed, main


@pytest.fixture
def mock_repeatmasker_content():
    """
    Mock RepeatMasker .out content based on the real sample1.fasta.out.
    First three lines are headers.
    """
    return (
        "   SW   perc perc perc  query           position in query           matching         repeat           position in repeat\n"
        "score   div. del. ins.  sequence        begin   end        (left)   repeat           class/family   begin  end    (left)   ID\n"
        "\n"
        "   19   12.4  0.0  2.8  NW_017385987.1     1405    1441 (4175035) + (TAT)n           Simple_repeat       1     36    (0)    1\n"
        " 8462    3.8  0.6  0.6  NW_017385987.1     1442    2500 (4173976) C rnd-1_family-30  LINE/Penelope   (260)   1059      1    2\n"
        " 7898    3.4  0.6  0.0  NW_017385987.1     3030    3991 (4172485) + rnd-1_family-30  LINE/Penelope     349   1316    (3)    3\n"
        "   14   22.9  0.0  6.1  NW_017385987.1     4643    4694 (4171782) + (TAGAGT)n        Simple_repeat       1     49    (0)    4\n"
    )


def test_parse_repeatmasker_out(mock_repeatmasker_content):
    mock_file = StringIO(mock_repeatmasker_content)
    df = parse_repeatmasker_out(mock_file)
    
    # Verify number of columns (15 as per the list in extract_te.py)
    assert len(df.columns) == 15
    
    # Verify number of rows (should skip 3 header lines)
    assert len(df) == 4
    
    # Verify specific content
    assert df.iloc[0]["query"] == "NW_017385987.1"
    assert df.iloc[0]["q_start"] == 1405
    assert df.iloc[1]["class_family"] == "LINE/Penelope"
    assert df.iloc[1]["strand"] == "C"


def test_filter_by_te_name(mock_repeatmasker_content):
    mock_file = StringIO(mock_repeatmasker_content)
    df = parse_repeatmasker_out(mock_file)
    
    # Filter for Penelope (case-insensitive)
    filtered = filter_by_te_name(df, "penelope")
    assert not filtered.empty
    assert len(filtered) == 2
    assert all(filtered["class_family"].str.contains("Penelope", case=False))
    
    # Filter for something non-existent
    filtered_empty = filter_by_te_name(df, "NON_EXISTENT")
    assert filtered_empty.empty


def test_convert_to_bed(mock_repeatmasker_content):
    mock_file = StringIO(mock_repeatmasker_content)
    df = parse_repeatmasker_out(mock_file)
    
    bed_df = convert_to_bed(df)
    
    # Verify BED columns: chrom, start, end, name, score, strand
    assert list(bed_df.columns) == ["chrom", "start", "end", "name", "score", "strand"]
    
    # Verify 1-based to 0-based conversion for the first row
    # RepeatMasker: 1405 -> BED: 1404
    assert bed_df.iloc[0]["start"] == 1404
    assert bed_df.iloc[0]["end"] == 1441
    
    # Verify strand conversion
    # Row 0: + -> +
    assert bed_df.iloc[0]["strand"] == "+"
    # Row 1: C -> -
    assert bed_df.iloc[1]["strand"] == "-"


def test_write_bed(tmp_path, mock_repeatmasker_content):
    mock_file = StringIO(mock_repeatmasker_content)
    df = parse_repeatmasker_out(mock_file)
    bed_df = convert_to_bed(df)
    
    output_bed = tmp_path / "test.bed"
    write_bed(bed_df, str(output_bed))
    
    assert output_bed.exists()
    
    # Read back and verify
    lines = output_bed.read_text().splitlines()
    assert len(lines) == 4
    # Check first line: chrom, start, end, name, score, strand
    # NW_017385987.1 1404 1441 repeat_1 19 +
    parts = lines[0].split("\t")
    assert parts[0] == "NW_017385987.1"
    assert parts[1] == "1404"
    assert parts[2] == "1441"
    assert parts[3] == "repeat_1"
    assert parts[4] == "19.0" or parts[4] == "19" # pandas to_csv might add .0
    assert parts[5] == "+"


def test_main_functionality(tmp_path, mock_repeatmasker_content):
    input_rm = tmp_path / "sample.fasta.out"
    output_bed = tmp_path / "output.bed"
    log_file = tmp_path / "test.log"
    
    input_rm.write_text(mock_repeatmasker_content)
    
    main(str(input_rm), str(output_bed), "Penelope", str(log_file))
    
    assert output_bed.exists()
    assert log_file.exists()
    
    lines = output_bed.read_text().splitlines()
    assert len(lines) == 2 # Penelope occurs twice in the mock data
