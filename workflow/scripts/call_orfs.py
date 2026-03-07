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


RANGE_PATTERN = re.compile(r"\[(?P<start>\d+)\s*-\s*(?P<end>\d+)\]")
REVERSE_PATTERN = re.compile(r"REVERSE SENSE", re.IGNORECASE)


def parse_getorf_header(record: SeqRecord):
    """Extract ORF metadata from a getorf (EMBOSS) FASTA record header.

    getorf appends ``_N`` (ORF number) to the sequence id and places coordinates
    as ``[start - end]`` followed by ``(REVERSE SENSE)`` for minus-strand ORFs.
    Missing fields are returned as ``None``.
    """
    header = record.description
    token = header.split()[0]
    # getorf appends _N to the original seqid; strip it to recover the TE id
    seqid = re.sub(r"_\d+$", "", token)
    start = None
    end = None

    m = RANGE_PATTERN.search(header)
    if m:
        start = int(m.group("start"))
        end = int(m.group("end"))

    strand = "-" if REVERSE_PATTERN.search(header) else "+"

    return {
        "te_id": seqid,
        "start": start,
        "end": end,
        "frame": None,
        "strand": strand,
    }


def run_getorf(input_fasta, out_path, find, min_orf_nt):
    """Run EMBOSS getorf on an input FASTA.

    Args:
        input_fasta: Path to nucleotide FASTA to scan.
        out_path: Path where getorf writes its output.
        find: getorf ``-find`` mode controlling which ORFs are reported:

            * ``0`` — protein translation of regions between stop codons.
              No start codon required; captures truncated (5'-incomplete) ORFs.
            * ``1`` — protein translation of canonical ORFs (start → stop).
              Requires a standard start codon (ATG).
            * ``2`` — nucleotide sequence of regions between stop codons.
              No start codon required; paired with ``find=0``.
            * ``3`` — nucleotide sequence of canonical ORFs (start → stop).
              Requires a standard start codon (ATG); paired with ``find=1``.

        min_orf_nt: Minimum ORF size in nucleotides.
    """
    cmd = [
        "getorf",
        "-sequence",
        input_fasta,
        "-outseq",
        out_path,
        "-minsize",
        str(min_orf_nt),
        "-find",
        str(find),
    ]
    logger.info("Running getorf: %s", " ".join(cmd))
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
    require_start_codon=True,
    log_file=None,
):
    """Call ORFs with getorf (EMBOSS) and export table + nucleotide/protein FASTAs.

    The function runs getorf twice (nt and aa outputs), pairs records by
    index, normalizes ORF ids, captures basic metadata for downstream workflow
    steps, and writes Snakemake outputs.

    Args:
        require_start_codon: If ``True``, only canonical ORFs starting with ATG
            are reported (getorf ``-find 1/3``).  If ``False``, any region
            between two stop codons is reported, capturing 5'-truncated ORFs
            (getorf ``-find 0/2``).
    """
    setup_logging(log_file, __name__)
    if require_start_codon:
        find_nt, find_aa = 3, 1
        logger.info("ORF mode: canonical (start→stop, ATG required)")
    else:
        find_nt, find_aa = 2, 0
        logger.info("ORF mode: permissive (stop→stop, any sense codon allowed)")

    min_orf_nt = int(min_orf_aa) * 3
    with tempfile.TemporaryDirectory(prefix="getorf_") as tmpdir:
        nt_tmp = str(Path(tmpdir) / "orfs_nt.fa")
        aa_tmp = str(Path(tmpdir) / "orfs_aa.fa")
        run_getorf(input_fasta, nt_tmp, find=find_nt, min_orf_nt=min_orf_nt)
        run_getorf(input_fasta, aa_tmp, find=find_aa, min_orf_nt=min_orf_nt)

        nt_records_raw = list(SeqIO.parse(nt_tmp, "fasta"))
        aa_records_raw = list(SeqIO.parse(aa_tmp, "fasta"))

    if len(nt_records_raw) != len(aa_records_raw):
        logger.warning(
            "getorf nucleotide/protein ORF count mismatch (%d vs %d). Pairing by minimum count.",
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
        parsed = parse_getorf_header(nt_rec)
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
                "has_start": require_start_codon,
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
            snakemake.params.require_start_codon,
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
        parser.add_argument(
            "--require-start-codon",
            action=argparse.BooleanOptionalAction,
            default=True,
            help="Require ATG start codon (default). Use --no-require-start-codon "
            "to include truncated ORFs starting at any sense codon.",
        )
        parser.add_argument("--log-file")
        args = parser.parse_args()
        main(
            args.input_fasta,
            args.output_table,
            args.output_nt_fasta,
            args.output_aa_fasta,
            args.min_orf_aa,
            args.require_start_codon,
            args.log_file,
        )
