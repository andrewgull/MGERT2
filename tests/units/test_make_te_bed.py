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

from workflow.scripts.make_te_bed import (
    parse_repeatmasker_out,
    filter_by_te_name,
    convert_to_bed,
    write_bed,
    find_out_file,
    main,
)


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
    
    # Verify number of columns (15 as per the list in make_te_bed.py)
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


def test_find_out_file_multiple(tmp_path):
    """
    Test find_out_file when multiple .out files are present.
    """
    rm_dir = tmp_path / "rm_outputs"
    rm_dir.mkdir()
    
    # Create multiple .out files
    file1 = rm_dir / "sample_A.fasta.out"
    file2 = rm_dir / "sample_B.fasta.out"
    file1.write_text("dummy A")
    file2.write_text("dummy B")
    
    # 1. No sample provided, should pick the first one (alphabetical usually)
    # The order of os.listdir can vary, but find_out_file takes files[0]
    result = find_out_file(str(rm_dir))
    assert result.endswith(".out")
    
    # 2. Sample provided, should pick the matching one
    result_a = find_out_file(str(rm_dir), sample="sample_A")
    assert os.path.basename(result_a) == "sample_A.fasta.out"
    
    result_b = find_out_file(str(rm_dir), sample="sample_B")
    assert os.path.basename(result_b) == "sample_B.fasta.out"
    
    # 3. Sample provided but no exact match for sample.fasta.out, 
    # but starts with sample
    file3 = rm_dir / "sample_C.out"
    file3.write_text("dummy C")
    result_c = find_out_file(str(rm_dir), sample="sample_C")
    assert os.path.basename(result_c) == "sample_C.out"


def test_find_out_file_single(tmp_path):
    rm_dir = tmp_path / "rm"
    rm_dir.mkdir()
    (rm_dir / "sample.fasta.out").write_text("dummy")
    result = find_out_file(str(rm_dir))
    assert result.endswith("sample.fasta.out")


def test_find_out_file_not_found(tmp_path):
    rm_dir = tmp_path / "rm"
    rm_dir.mkdir()
    with pytest.raises(FileNotFoundError):
        find_out_file(str(rm_dir))


def test_parse_repeatmasker_out_bad_line(tmp_path, mock_repeatmasker_content, caplog):
    # A line with 16 tokens (one too many for the 15-column schema) triggers on_bad_line
    bad_line = "   19   12.4  0.0  2.8  NW_017385987.1     1405    1441 (4175035) + (TAT)n           Simple_repeat       1     36    (0)    1    EXTRA\n"
    bad_file = tmp_path / "bad.out"
    bad_file.write_text(mock_repeatmasker_content + bad_line)
    with caplog.at_level("WARNING"):
        df = parse_repeatmasker_out(str(bad_file))
    assert "Skipping bad line" in caplog.text


def test_main_no_hits_logs_warning(tmp_path, mock_repeatmasker_content, caplog):
    out_file = tmp_path / "sample.fasta.out"
    out_file.write_text(mock_repeatmasker_content)
    bed_out = tmp_path / "out.bed"
    with caplog.at_level("WARNING"):
        main(str(out_file), str(bed_out), "NONEXISTENT_TE", None)
    assert "No hits found" in caplog.text


def test_convert_to_bed_unexpected_strand(caplog):
    """
    Test convert_to_bed when unexpected strand values are present.
    It should log a warning.
    """
    data = {
        "query": ["seq1"],
        "q_start": [100],
        "q_end": [200],
        "id": [1],
        "sw_score": [1000],
        "strand": ["?"],  # Unexpected strand
    }
    df = pd.DataFrame(data)
    
    with caplog.at_level("WARNING"):
        bed_df = convert_to_bed(df)
        
    assert "Found 1 rows with unexpected strand values" in caplog.text
    assert bed_df.iloc[0]["strand"] == "?"
