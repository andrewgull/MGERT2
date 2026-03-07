import logging
import subprocess
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from utils import setup_logging

logger = logging.getLogger(__name__)

RPS_COLS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "stitle",
]


def write_empty_outputs(aa_fasta, hits_out, summary_out):
    orf_ids = [rec.id for rec in SeqIO.parse(aa_fasta, "fasta")]
    pd.DataFrame(columns=RPS_COLS + ["qcov"]).to_csv(hits_out, sep="\t", index=False)
    summary = pd.DataFrame(
        {
            "orf_id": orf_ids,
            "domain_support": [False] * len(orf_ids),
            "hit_count": [0] * len(orf_ids),
            "best_domain": [""] * len(orf_ids),
            "best_evalue": [float("nan")] * len(orf_ids),
            "best_bitscore": [float("nan")] * len(orf_ids),
        }
    )
    summary.to_csv(summary_out, sep="\t", index=False)


def main(
    aa_fasta,
    hits_out,
    summary_out,
    db,
    evalue,
    max_target_seqs,
    min_query_coverage,
    min_bitscore,
    threads=1,
    log_file=None,
):
    setup_logging(log_file, __name__)
    Path(hits_out).parent.mkdir(parents=True, exist_ok=True)

    if not db:
        logger.warning("No rpsblast database configured; writing empty domain outputs")
        write_empty_outputs(aa_fasta, hits_out, summary_out)
        return

    db_candidates = list(Path().glob(f"{db}*")) if not Path(db).exists() else [Path(db)]
    if not db_candidates:
        logger.warning(
            "Configured rpsblast database '%s' not found; writing empty domain outputs",
            db,
        )
        write_empty_outputs(aa_fasta, hits_out, summary_out)
        return

    raw_out = str(Path(hits_out).with_suffix(".raw.tsv"))
    cmd = [
        "rpsblast",
        "-query",
        aa_fasta,
        "-db",
        db,
        "-evalue",
        str(evalue),
        "-max_target_seqs",
        str(max_target_seqs),
        "-num_threads",
        str(threads),
        "-outfmt",
        "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle",
        "-out",
        raw_out,
    ]

    logger.info("Running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)

    qlen = {rec.id: len(rec.seq) for rec in SeqIO.parse(aa_fasta, "fasta")}

    if Path(raw_out).stat().st_size == 0:
        logger.info("No RPS-BLAST hits found")
        write_empty_outputs(aa_fasta, hits_out, summary_out)
        return

    hits = pd.read_csv(raw_out, sep="\t", names=RPS_COLS)
    hits["qcov"] = hits.apply(lambda r: float(r["length"]) / max(qlen.get(r["qseqid"], 1), 1), axis=1)
    hits = hits[
        (hits["evalue"] <= float(evalue))
        & (hits["qcov"] >= float(min_query_coverage))
        & (hits["bitscore"] >= float(min_bitscore))
    ]
    hits.to_csv(hits_out, sep="\t", index=False)

    if hits.empty:
        write_empty_outputs(aa_fasta, hits_out, summary_out)
        return

    best = (
        hits.sort_values(["qseqid", "bitscore", "evalue"], ascending=[True, False, True])
        .groupby("qseqid", as_index=False)
        .first()
    )
    counts = hits.groupby("qseqid").size().rename("hit_count").reset_index()
    summary = best.merge(counts, on="qseqid", how="left")
    summary = summary.rename(
        columns={
            "qseqid": "orf_id",
            "stitle": "best_domain",
            "evalue": "best_evalue",
            "bitscore": "best_bitscore",
        }
    )[["orf_id", "hit_count", "best_domain", "best_evalue", "best_bitscore"]]
    summary["domain_support"] = True

    all_orfs = pd.DataFrame({"orf_id": list(qlen.keys())})
    summary = all_orfs.merge(summary, on="orf_id", how="left")
    summary["domain_support"] = summary["domain_support"].fillna(False)
    summary["hit_count"] = summary["hit_count"].fillna(0).astype(int)
    summary["best_domain"] = summary["best_domain"].fillna("")
    summary.to_csv(summary_out, sep="\t", index=False)

    logger.info("Wrote %d filtered hits and %d ORF domain summaries", len(hits), len(summary))


if __name__ == "__main__":
    try:
        main(
            snakemake.input.aa_fasta,
            snakemake.output.hits,
            snakemake.output.summary,
            snakemake.params.db,
            snakemake.params.evalue,
            snakemake.params.max_target_seqs,
            snakemake.params.min_query_coverage,
            snakemake.params.min_bitscore,
            snakemake.threads,
            snakemake.log[0] if snakemake.log else None,
        )
    except NameError:
        import argparse

        parser = argparse.ArgumentParser(description="Run RPS-BLAST on ORFs")
        parser.add_argument("--aa-fasta", required=True)
        parser.add_argument("--hits-out", required=True)
        parser.add_argument("--summary-out", required=True)
        parser.add_argument("--db", default="")
        parser.add_argument("--evalue", type=float, default=1e-5)
        parser.add_argument("--max-target-seqs", type=int, default=10)
        parser.add_argument("--min-query-coverage", type=float, default=0.35)
        parser.add_argument("--min-bitscore", type=float, default=50.0)
        parser.add_argument("--threads", type=int, default=1)
        parser.add_argument("--log-file")
        args = parser.parse_args()
        main(
            args.aa_fasta,
            args.hits_out,
            args.summary_out,
            args.db,
            args.evalue,
            args.max_target_seqs,
            args.min_query_coverage,
            args.min_bitscore,
            args.threads,
            args.log_file,
        )
