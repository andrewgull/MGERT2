import logging
from pathlib import Path

import pandas as pd
from utils import setup_logging

logger = logging.getLogger(__name__)


def main(raw_orfs_tsv, classified_tsv, summary_tsv, summary_html, log_file=None):
    setup_logging(log_file, __name__)

    raw_df = pd.read_csv(raw_orfs_tsv, sep="\t")
    cls_df = pd.read_csv(classified_tsv, sep="\t")

    total_raw_orfs = len(raw_df)
    total_filtered_orfs = len(cls_df)
    te_with_raw_orfs = raw_df["te_id"].nunique() if "te_id" in raw_df.columns else 0
    te_with_filtered_orfs = (
        cls_df["te_id"].nunique() if "te_id" in cls_df.columns else 0
    )
    longest_orf_aa = int(cls_df["aa_len"].max()) if not cls_df.empty else 0

    class_counts = (
        cls_df["confidence_class"]
        .value_counts(dropna=False)
        .rename_axis("metric")
        .reset_index(name="value")
        if "confidence_class" in cls_df.columns
        else pd.DataFrame(columns=["metric", "value"])
    )

    core_metrics = pd.DataFrame(
        [
            {"metric": "total_raw_orfs", "value": total_raw_orfs},
            {"metric": "total_filtered_orfs", "value": total_filtered_orfs},
            {"metric": "te_with_raw_orfs", "value": te_with_raw_orfs},
            {"metric": "te_with_filtered_orfs", "value": te_with_filtered_orfs},
            {"metric": "longest_orf_aa", "value": longest_orf_aa},
        ]
    )

    out = pd.concat([core_metrics, class_counts], ignore_index=True)
    Path(summary_tsv).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(summary_tsv, sep="\t", index=False)

    top_orfs = (
        cls_df.sort_values("aa_len", ascending=False).head(25)
        if "aa_len" in cls_df.columns
        else cls_df
    )
    html = ["<html><body>", "<h1>ORF Coding Potential Report</h1>"]
    html.append("<h2>Summary Metrics</h2>")
    html.append(out.to_html(index=False))
    html.append("<h2>Top ORFs</h2>")
    html.append(top_orfs.to_html(index=False))
    html.append("</body></html>")
    Path(summary_html).write_text("\n".join(html))

    logger.info("Wrote ORF summary report to %s and %s", summary_tsv, summary_html)


if __name__ == "__main__":
    try:
        main(
            snakemake.input.raw_orfs,
            snakemake.input.classified,
            snakemake.output.summary_tsv,
            snakemake.output.summary_html,
            snakemake.log[0] if snakemake.log else None,
        )
    except NameError:
        import argparse

        parser = argparse.ArgumentParser(
            description="Create ORF coding potential report"
        )
        parser.add_argument("--raw-orfs", required=True)
        parser.add_argument("--classified", required=True)
        parser.add_argument("--summary-tsv", required=True)
        parser.add_argument("--summary-html", required=True)
        parser.add_argument("--log-file")
        args = parser.parse_args()
        main(
            args.raw_orfs,
            args.classified,
            args.summary_tsv,
            args.summary_html,
            args.log_file,
        )
