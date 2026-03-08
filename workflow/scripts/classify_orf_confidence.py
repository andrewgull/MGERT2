import logging
from pathlib import Path

import pandas as pd
from utils import setup_logging

logger = logging.getLogger(__name__)


def classify_row(r, high_min_aa, high_min_intrinsic, putative_min_aa, putative_min_intrinsic):
    aa_len = float(r.get("aa_len", 0))
    intrinsic = float(r.get("intrinsic_score", 0))
    domain_support = bool(r.get("domain_support", False))

    if aa_len >= high_min_aa and intrinsic >= high_min_intrinsic and domain_support:
        return "high_confidence_coding"
    if aa_len >= putative_min_aa and (intrinsic >= putative_min_intrinsic or domain_support):
        return "putative_coding"
    return "unlikely_coding"


def main(
    filtered_orfs,
    intrinsic_scores,
    domain_summary,
    output_tsv,
    high_min_aa,
    high_min_intrinsic,
    putative_min_aa,
    putative_min_intrinsic,
    log_file=None,
):
    setup_logging(log_file, __name__)

    orf_df = pd.read_csv(filtered_orfs, sep="\t")
    intrinsic_df = pd.read_csv(intrinsic_scores, sep="\t")
    domain_df = pd.read_csv(domain_summary, sep="\t")

    merged = orf_df.merge(intrinsic_df, on="orf_id", how="left").merge(domain_df, on="orf_id", how="left")
    merged["domain_support"] = merged["domain_support"].fillna(False)
    merged["intrinsic_score"] = merged["intrinsic_score"].fillna(0.0) if "intrinsic_score" in merged.columns else 0.0
    merged["aa_len"] = merged["aa_len"].fillna(0.0) if "aa_len" in merged.columns else 0.0
    merged["confidence_class"] = merged.apply(
        classify_row,
        axis=1,
        args=(high_min_aa, high_min_intrinsic, putative_min_aa, putative_min_intrinsic),
    )

    cols = [
        "te_id",
        "orf_id",
        "aa_len",
        "intrinsic_score",
        "intrinsic_label",
        "domain_support",
        "hit_count",
        "best_domain",
        "best_evalue",
        "best_bitscore",
        "confidence_class",
    ]
    for c in cols:
        if c not in merged.columns:
            merged[c] = pd.NA

    Path(output_tsv).parent.mkdir(parents=True, exist_ok=True)
    merged[cols].to_csv(output_tsv, sep="\t", index=False)
    logger.info("Classified %d ORFs", len(merged))


if __name__ == "__main__":
    if "snakemake" in globals():
        main(
            snakemake.input.filtered_orfs,
            snakemake.input.intrinsic_scores,
            snakemake.input.domain_summary,
            snakemake.output.tsv,
            snakemake.params.high_min_aa,
            snakemake.params.high_min_intrinsic,
            snakemake.params.putative_min_aa,
            snakemake.params.putative_min_intrinsic,
            snakemake.log[0] if snakemake.log else None,
        )
    else:
        import argparse

        parser = argparse.ArgumentParser(description="Classify ORF coding confidence")
        parser.add_argument("--filtered-orfs", required=True)
        parser.add_argument("--intrinsic-scores", required=True)
        parser.add_argument("--domain-summary", required=True)
        parser.add_argument("--output-tsv", required=True)
        parser.add_argument("--high-min-aa", type=float, default=300)
        parser.add_argument("--high-min-intrinsic", type=float, default=0.6)
        parser.add_argument("--putative-min-aa", type=float, default=150)
        parser.add_argument("--putative-min-intrinsic", type=float, default=0.5)
        parser.add_argument("--log-file")
        args = parser.parse_args()
        main(
            args.filtered_orfs,
            args.intrinsic_scores,
            args.domain_summary,
            args.output_tsv,
            args.high_min_aa,
            args.high_min_intrinsic,
            args.putative_min_aa,
            args.putative_min_intrinsic,
            args.log_file,
        )
