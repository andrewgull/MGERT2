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


def load_domain_map(domains_csv):
    """Load mapping from .smp stem to domain type from a two-column TSV (filename.smp, domain_type).

    Returns an empty dict if the file is missing or not provided.
    """
    if not domains_csv or not Path(domains_csv).exists():
        return {}
    df = pd.read_csv(
        domains_csv, sep="\t", header=None, names=["filename", "domain_type"]
    )
    return {Path(row.filename).stem: row.domain_type for _, row in df.iterrows()}


def get_domain_type(stitle, domain_map):
    """Extract domain type from an RPS-BLAST stitle using the domain map.

    stitle is typically "pfam00078.1, RVT_1, Reverse transcriptase...".
    The accession is the first comma-separated token with any version suffix stripped.
    """
    if not domain_map or not isinstance(stitle, str):
        return ""
    accession = stitle.split(",")[0].strip()
    if "." in accession:
        accession = accession.rsplit(".", 1)[0]
    return domain_map.get(accession, "")


def _domain_types_for_group(series):
    """Return pipe-separated sorted unique non-empty domain types from a Series."""
    return "|".join(sorted(set(t for t in series if t)))


def _n_domain_types_for_group(series):
    return len(set(t for t in series if t))


def write_empty_outputs(aa_fasta, hits_out, summary_out):
    orf_ids = [rec.id for rec in SeqIO.parse(aa_fasta, "fasta")]
    pd.DataFrame(columns=RPS_COLS + ["qcov", "domain_type"]).to_csv(
        hits_out, sep="\t", index=False
    )
    summary = pd.DataFrame(
        {
            "orf_id": orf_ids,
            "domain_support": [False] * len(orf_ids),
            "hit_count": [0] * len(orf_ids),
            "best_domain": [""] * len(orf_ids),
            "best_evalue": [float("nan")] * len(orf_ids),
            "best_bitscore": [float("nan")] * len(orf_ids),
            "domain_types": [""] * len(orf_ids),
            "n_domain_types": [0] * len(orf_ids),
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
    domains_csv="",
    threads=1,
    log_file=None,
):
    setup_logging(log_file, __name__)
    Path(hits_out).parent.mkdir(parents=True, exist_ok=True)

    domain_map = load_domain_map(domains_csv)
    if domain_map:
        logger.info(
            "Loaded %d domain type mappings from %s", len(domain_map), domains_csv
        )
    else:
        logger.info("No domain type map loaded; domain_types column will be empty")

    if not db:
        logger.warning("No rpsblast database configured; writing empty domain outputs")
        write_empty_outputs(aa_fasta, hits_out, summary_out)
        return

    db_path = Path(db)
    if db_path.exists():
        db_candidates = [db_path]
    elif db_path.is_absolute():
        db_candidates = list(db_path.parent.glob(f"{db_path.name}*"))
    else:
        db_candidates = list(Path().glob(f"{db}*"))
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
    hits["qcov"] = hits.apply(
        lambda r: float(r["qend"] - r["qstart"] + 1) / max(qlen.get(r["qseqid"], 1), 1),
        axis=1,
    )
    hits = hits[
        (hits["evalue"] <= float(evalue))
        & (hits["qcov"] >= float(min_query_coverage))
        & (hits["bitscore"] >= float(min_bitscore))
    ]
    hits["domain_type"] = hits["stitle"].apply(lambda s: get_domain_type(s, domain_map))
    hits.to_csv(hits_out, sep="\t", index=False)

    if hits.empty:
        write_empty_outputs(aa_fasta, hits_out, summary_out)
        return

    best = (
        hits.sort_values(
            ["qseqid", "bitscore", "evalue"], ascending=[True, False, True]
        )
        .groupby("qseqid", as_index=False)
        .first()
    )
    counts = hits.groupby("qseqid").size().rename("hit_count").reset_index()
    type_agg = (
        hits.groupby("qseqid")["domain_type"]
        .apply(_domain_types_for_group)
        .rename("domain_types")
        .reset_index()
    )
    n_type_agg = (
        hits.groupby("qseqid")["domain_type"]
        .apply(_n_domain_types_for_group)
        .rename("n_domain_types")
        .reset_index()
    )

    summary = best.merge(counts, on="qseqid", how="left")
    summary = summary.merge(type_agg, on="qseqid", how="left")
    summary = summary.merge(n_type_agg, on="qseqid", how="left")
    summary = summary.rename(
        columns={
            "qseqid": "orf_id",
            "stitle": "best_domain",
            "evalue": "best_evalue",
            "bitscore": "best_bitscore",
        }
    )[
        [
            "orf_id",
            "hit_count",
            "best_domain",
            "best_evalue",
            "best_bitscore",
            "domain_types",
            "n_domain_types",
        ]
    ]
    summary["domain_support"] = True

    all_orfs = pd.DataFrame({"orf_id": list(qlen.keys())})
    summary = all_orfs.merge(summary, on="orf_id", how="left")
    summary["domain_support"] = summary["domain_support"].fillna(False)
    summary["hit_count"] = summary["hit_count"].fillna(0).astype(int)
    summary["best_domain"] = summary["best_domain"].fillna("")
    summary["domain_types"] = summary["domain_types"].fillna("")
    summary["n_domain_types"] = summary["n_domain_types"].fillna(0).astype(int)
    summary.to_csv(summary_out, sep="\t", index=False)

    logger.info(
        "Wrote %d filtered hits and %d ORF domain summaries", len(hits), len(summary)
    )


if __name__ == "__main__":
    if "snakemake" in globals():
        main(
            snakemake.input.aa_fasta,
            snakemake.output.hits,
            snakemake.output.summary,
            snakemake.params.db,
            snakemake.params.evalue,
            snakemake.params.max_target_seqs,
            snakemake.params.min_query_coverage,
            snakemake.params.min_bitscore,
            snakemake.params.domains_csv,
            snakemake.threads,
            snakemake.log[0] if snakemake.log else None,
        )
    else:
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
        parser.add_argument("--domains-csv", default="")
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
            args.domains_csv,
            args.threads,
            args.log_file,
        )
