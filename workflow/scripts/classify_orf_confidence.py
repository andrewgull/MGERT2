import logging
from pathlib import Path

import pandas as pd
from utils import setup_logging

logger = logging.getLogger(__name__)


def classify_row(
    r, high_min_aa, high_min_intrinsic, putative_min_aa, putative_min_intrinsic
):
    aa_len = float(r.get("aa_len", 0))
    intrinsic = float(r.get("intrinsic_score", 0))
    domain_support = bool(r.get("domain_support", False))

    if aa_len >= high_min_aa and intrinsic >= high_min_intrinsic and domain_support:
        return "high_confidence_coding"
    if aa_len >= putative_min_aa and (
        intrinsic >= putative_min_intrinsic or domain_support
    ):
        return "putative_coding"
    return "unlikely_coding"


def _parse_domain_types(s):
    """Return the set of domain types from a pipe-separated string."""
    if not isinstance(s, str) or not s:
        return set()
    return set(s.split("|"))


def apply_domain_filter(
    merged,
    required_domains,
    min_domain_types_per_orf,
    allow_split_across_orfs,
    te_coverage_scope,
):
    """Add per-ORF and per-TE domain filter columns to *merged* (in-place) and return it.

    Columns added
    -------------
    passes_domain_filter : bool
        True when the ORF has at least *min_domain_types_per_orf* distinct types
        that appear in *required_domains* (or any types when *required_domains* is empty).
    te_domain_complete : bool
        When *allow_split_across_orfs* is False: True when the TE has at least one
        ORF with ``passes_domain_filter=True``.
        When *allow_split_across_orfs* is True: True when the union of domain_types
        across the TE's ORFs covers every entry in *required_domains*.
        The set of ORFs used for the union is controlled by *te_coverage_scope*:
        ``"all"`` includes every ORF; ``"passing"`` includes only ORFs that already
        pass the per-ORF filter.
        When *required_domains* is empty, ``te_domain_complete`` is always True.
    """
    req = set(required_domains) if required_domains else set()

    # --- per-ORF filter ---
    def _per_orf_passes(domain_types_str):
        types = _parse_domain_types(domain_types_str)
        if not req:
            return len(types) >= min_domain_types_per_orf
        return len(types & req) >= min_domain_types_per_orf

    domain_types_col = (
        merged["domain_types"]
        if "domain_types" in merged.columns
        else pd.Series("", index=merged.index)
    )
    merged["passes_domain_filter"] = domain_types_col.apply(_per_orf_passes)

    # --- per-TE coverage ---
    if not req:
        merged["te_domain_complete"] = True
    elif allow_split_across_orfs:
        if te_coverage_scope == "passing":
            source = merged[merged["passes_domain_filter"]]
        else:  # "all"
            source = merged
        te_coverage = source.groupby("te_id")["domain_types"].apply(
            lambda series: set().union(*(_parse_domain_types(s) for s in series))
        )
        te_complete = te_coverage.apply(lambda types: req <= types)
        merged["te_domain_complete"] = merged["te_id"].map(te_complete).fillna(False)
    else:
        # Non-split mode: TE passes when at least one ORF satisfies the per-ORF filter.
        te_any_pass = merged.groupby("te_id")["passes_domain_filter"].any()
        merged["te_domain_complete"] = merged["te_id"].map(te_any_pass).fillna(False)

    return merged


def main(
    filtered_orfs,
    intrinsic_scores,
    domain_summary,
    output_tsv,
    high_min_aa,
    high_min_intrinsic,
    putative_min_aa,
    putative_min_intrinsic,
    required_domains=None,
    min_domain_types_per_orf=1,
    allow_split_across_orfs=False,
    te_coverage_scope="all",
    log_file=None,
):
    setup_logging(log_file, __name__)

    orf_df = pd.read_csv(filtered_orfs, sep="\t")
    intrinsic_df = pd.read_csv(intrinsic_scores, sep="\t")
    domain_df = pd.read_csv(domain_summary, sep="\t")

    merged = orf_df.merge(intrinsic_df, on="orf_id", how="left").merge(
        domain_df, on="orf_id", how="left"
    )
    merged["domain_support"] = merged["domain_support"].fillna(False)
    merged["intrinsic_score"] = (
        merged["intrinsic_score"].fillna(0.0)
        if "intrinsic_score" in merged.columns
        else 0.0
    )
    merged["aa_len"] = (
        merged["aa_len"].fillna(0.0) if "aa_len" in merged.columns else 0.0
    )
    merged["domain_types"] = (
        merged["domain_types"].fillna("") if "domain_types" in merged.columns else ""
    )
    merged["n_domain_types"] = (
        merged["n_domain_types"].fillna(0) if "n_domain_types" in merged.columns else 0
    )

    merged["confidence_class"] = merged.apply(
        classify_row,
        axis=1,
        args=(high_min_aa, high_min_intrinsic, putative_min_aa, putative_min_intrinsic),
    )

    merged = apply_domain_filter(
        merged,
        required_domains or [],
        min_domain_types_per_orf,
        allow_split_across_orfs,
        te_coverage_scope,
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
        "domain_types",
        "n_domain_types",
        "passes_domain_filter",
        "te_domain_complete",
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
            snakemake.params.required_domains,
            snakemake.params.min_domain_types_per_orf,
            snakemake.params.allow_split_across_orfs,
            snakemake.params.te_coverage_scope,
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
        parser.add_argument("--required-domains", nargs="*", default=[])
        parser.add_argument("--min-domain-types-per-orf", type=int, default=1)
        parser.add_argument("--allow-split-across-orfs", action="store_true")
        parser.add_argument(
            "--te-coverage-scope", choices=["all", "passing"], default="all"
        )
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
            args.required_domains,
            args.min_domain_types_per_orf,
            args.allow_split_across_orfs,
            args.te_coverage_scope,
            args.log_file,
        )
