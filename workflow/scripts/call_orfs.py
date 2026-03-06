import logging
import re
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from utils import setup_logging

logger = logging.getLogger(__name__)


RANGE_PATTERNS = [
    re.compile(r"(?P<seqid>\S+):(?P<start>\d+)-(?P<end>\d+)"),
    re.compile(r"\[(?P<start>\d+)\s*-\s*(?P<end>\d+)\]"),
]
FRAME_PATTERN = re.compile(r"frame(?:=|\s)(?P<frame>-?\d+)", re.IGNORECASE)


def parse_orffinder_header(record: SeqRecord):
    """Extract ORF metadata from an ORFfinder FASTA record header.

    The parser pulls sequence id, start/end coordinates, reading frame, and
    strand from ORFfinder's description string using permissive regex patterns.
    Missing fields are returned as ``None``.
    """
    header = record.description
    token = header.split()[0]
    seqid = token.split(":")[0]
    start = None
    end = None
    frame = None

    for pat in RANGE_PATTERNS:
        m = pat.search(header)
        if m:
            start = int(m.group("start"))
            end = int(m.group("end"))
            seqid = m.groupdict().get("seqid", seqid)
            break

    fm = FRAME_PATTERN.search(header)
    if fm:
        frame = int(fm.group("frame"))

    strand = "+"
    if frame is not None and frame < 0:
        strand = "-"

    return {
        "te_id": seqid,
        "start": start,
        "end": end,
        "frame": frame,
        "strand": strand,
    }


def run_orffinder(input_fasta, out_path, outfmt, min_orf_aa):
    """Run NCBI ORFfinder on an input FASTA with minimum ORF length in aa.

    Args:
        input_fasta: Path to nucleotide FASTA to scan.
        out_path: Path where ORFfinder writes its output.
        outfmt: ORFfinder output format integer (1=nt FASTA, 2=aa FASTA).
        min_orf_aa: Minimum ORF size in amino acids.
    """
    min_orf_nt = int(min_orf_aa) * 3
    cmd = [
        "ORFfinder",
        "-in",
        input_fasta,
        "-out",
        out_path,
        "-outfmt",
        str(outfmt),
        "-ml",
        str(min_orf_nt),
        "-strand",
        "both",
    ]
    logger.info("Running ORFfinder: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


def write_fasta(records, output_path):
    """Write SeqRecord objects to FASTA, creating parent directories if needed."""
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, output_path, "fasta")


def main(
    input_fasta,
    output_table,
    output_nt_fasta,
    output_aa_fasta,
    min_orf_aa,
    start_codons,
    stop_codons,
    log_file=None,
):
    """Call ORFs with ORFfinder and export table + nucleotide/protein FASTAs.

    The function runs ORFfinder twice (nt and aa outputs), pairs records by
    index, normalizes ORF ids, captures basic metadata for downstream workflow
    steps, and writes Snakemake outputs.
    """
    setup_logging(log_file, __name__)
    if start_codons != ["ATG"] or sorted(stop_codons) != ["TAA", "TAG", "TGA"]:
        logger.warning(
            "Configured start/stop codons are not directly applied by this ORFfinder wrapper. "
            "Using ORFfinder defaults for start/stop handling."
        )

    with tempfile.TemporaryDirectory(prefix="orffinder_") as tmpdir:
        nt_tmp = str(Path(tmpdir) / "orfs_nt.fa")
        aa_tmp = str(Path(tmpdir) / "orfs_aa.fa")
        run_orffinder(input_fasta, nt_tmp, outfmt=1, min_orf_aa=min_orf_aa)
        run_orffinder(input_fasta, aa_tmp, outfmt=2, min_orf_aa=min_orf_aa)

        nt_records_raw = list(SeqIO.parse(nt_tmp, "fasta"))
        aa_records_raw = list(SeqIO.parse(aa_tmp, "fasta"))

    if len(nt_records_raw) != len(aa_records_raw):
        logger.warning(
            "ORFfinder nucleotide/protein ORF count mismatch (%d vs %d). Pairing by minimum count.",
            len(nt_records_raw),
            len(aa_records_raw),
        )
    n = min(len(nt_records_raw), len(aa_records_raw))

    rows = []
    nt_records = []
    aa_records = []
    for idx in range(n):
        nt_rec = nt_records_raw[idx]
        aa_rec = aa_records_raw[idx]
        parsed = parse_orffinder_header(nt_rec)
        orf_id = f"{parsed['te_id']}|orf_{idx+1}"
        rows.append(
            {
                "te_id": parsed["te_id"],
                "orf_id": orf_id,
                "start": parsed["start"],
                "end": parsed["end"],
                "strand": parsed["strand"],
                "frame": parsed["frame"],
                "nt_len": len(nt_rec.seq),
                "aa_len": len(aa_rec.seq),
                "has_start": True,
                "has_stop": True,
            }
        )
        nt_rec.id = orf_id
        nt_rec.description = ""
        aa_rec.id = orf_id
        aa_rec.description = ""
        nt_records.append(nt_rec)
        aa_records.append(aa_rec)

    Path(output_table).parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    if df.empty:
        df = pd.DataFrame(
            columns=[
                "te_id",
                "orf_id",
                "start",
                "end",
                "strand",
                "frame",
                "nt_len",
                "aa_len",
                "has_start",
                "has_stop",
            ]
        )
    df.to_csv(output_table, sep="\t", index=False)
    write_fasta(nt_records, output_nt_fasta)
    write_fasta(aa_records, output_aa_fasta)
    logger.info("Detected %d ORFs across %s", len(rows), input_fasta)


if __name__ == "__main__":
    try:
        main(
            snakemake.input.fasta,
            snakemake.output.table,
            snakemake.output.nt_fasta,
            snakemake.output.aa_fasta,
            snakemake.params.min_orf_aa,
            snakemake.params.start_codons,
            snakemake.params.stop_codons,
            snakemake.log[0] if snakemake.log else None,
        )
    except NameError:
        import argparse

        parser = argparse.ArgumentParser(description="Call ORFs from TE fasta")
        parser.add_argument("--input-fasta", required=True)
        parser.add_argument("--output-table", required=True)
        parser.add_argument("--output-nt-fasta", required=True)
        parser.add_argument("--output-aa-fasta", required=True)
        parser.add_argument("--min-orf-aa", type=int, default=100)
        parser.add_argument("--start-codons", nargs="+", default=["ATG"])
        parser.add_argument("--stop-codons", nargs="+", default=["TAA", "TAG", "TGA"])
        parser.add_argument("--log-file")
        args = parser.parse_args()
        main(
            args.input_fasta,
            args.output_table,
            args.output_nt_fasta,
            args.output_aa_fasta,
            args.min_orf_aa,
            args.start_codons,
            args.stop_codons,
            args.log_file,
        )
