import logging
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from utils import setup_logging

logger = logging.getLogger(__name__)


def gc_content(seq):
    seq = seq.upper()
    if not seq:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq)


def codon_diversity(seq):
    codons = [seq[i:i + 3] for i in range(0, len(seq) - 2, 3)]
    codons = [c for c in codons if len(c) == 3]
    if not codons:
        return 0.0
    return len(set(codons)) / len(codons)


def intrinsic_label(score, high_threshold, medium_threshold):
    if score >= high_threshold:
        return "high"
    if score >= medium_threshold:
        return "medium"
    return "low"


def main(
    input_table,
    input_nt_fasta,
    output_tsv,
    high_threshold,
    medium_threshold,
    log_file=None,
):
    setup_logging(log_file, __name__)

    if not (0.0 <= medium_threshold <= 1.0 and 0.0 <= high_threshold <= 1.0):
        raise ValueError(f"Thresholds must be in [0, 1]: high={high_threshold}, medium={medium_threshold}")
    if medium_threshold > high_threshold:
        raise ValueError(f"medium_threshold ({medium_threshold}) must be <= high_threshold ({high_threshold})")

    df = pd.read_csv(input_table, sep="\t")
    seqs = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(input_nt_fasta, "fasta")}

    missing = set(df["orf_id"]) - set(seqs)
    if missing:
        raise ValueError(f"ORF IDs present in table but missing from FASTA: {sorted(missing)}")

    rows = []
    for _, row in df.iterrows():
        orf_id = row["orf_id"]
        nt = seqs[orf_id]
        aa_len = float(row.get("aa_len", 0))

        length_score = min(aa_len / 300.0, 1.0)
        gc = gc_content(nt)
        gc_score = max(0.0, 1.0 - (abs(gc - 0.5) / 0.5))
        cd = codon_diversity(nt)

        score = (0.5 * length_score) + (0.3 * cd) + (0.2 * gc_score)

        rows.append(
            {
                "orf_id": orf_id,
                "gc_content": round(gc, 4),
                "codon_diversity": round(cd, 4),
                "length_score": round(length_score, 4),
                "intrinsic_score": round(score, 4),
                "intrinsic_label": intrinsic_label(score, high_threshold, medium_threshold),
            }
        )

    out = pd.DataFrame(rows)
    if out.empty:
        out = pd.DataFrame(
            columns=[
                "orf_id",
                "gc_content",
                "codon_diversity",
                "length_score",
                "intrinsic_score",
                "intrinsic_label",
            ]
        )
    Path(output_tsv).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_tsv, sep="\t", index=False)
    logger.info("Wrote intrinsic coding potential for %d ORFs", len(out))


if __name__ == "__main__":
    if "snakemake" in globals():
        main(
            snakemake.input.table,
            snakemake.input.nt_fasta,
            snakemake.output.tsv,
            snakemake.params.high_threshold,
            snakemake.params.medium_threshold,
            snakemake.log[0] if snakemake.log else None,
        )
    else:
        import argparse

        parser = argparse.ArgumentParser(description="Compute intrinsic coding potential scores")
        parser.add_argument("--input-table", required=True)
        parser.add_argument("--input-nt-fasta", required=True)
        parser.add_argument("--output-tsv", required=True)
        parser.add_argument("--high-threshold", type=float, default=0.7)
        parser.add_argument("--medium-threshold", type=float, default=0.45)
        parser.add_argument("--log-file")
        args = parser.parse_args()
        main(
            args.input_table,
            args.input_nt_fasta,
            args.output_tsv,
            args.high_threshold,
            args.medium_threshold,
            args.log_file,
        )
