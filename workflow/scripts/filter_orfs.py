import logging
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from utils import setup_logging

logger = logging.getLogger(__name__)


def load_fasta_dict(path):
    return {rec.id: rec for rec in SeqIO.parse(path, "fasta")}


def write_subset_fasta(records, keep_ids, output_path):
    keep = [records[k] for k in keep_ids if k in records]
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(keep, output_path, "fasta")


def main(
    input_table,
    input_nt_fasta,
    input_aa_fasta,
    output_table,
    output_nt_fasta,
    output_aa_fasta,
    min_orf_aa,
    max_orfs_per_te,
    require_stop_codon=True,
    log_file=None,
):
    setup_logging(log_file, __name__)

    Path(output_table).parent.mkdir(parents=True, exist_ok=True)
    Path(output_nt_fasta).parent.mkdir(parents=True, exist_ok=True)
    Path(output_aa_fasta).parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_table, sep="\t")
    if df.empty:
        df.to_csv(output_table, sep="\t", index=False)
        Path(output_nt_fasta).write_text("")
        Path(output_aa_fasta).write_text("")
        logger.warning("No ORFs to filter; produced empty outputs")
        return

    filtered = df[df["aa_len"] >= int(min_orf_aa)].copy()
    if require_stop_codon and "has_stop" in filtered.columns:
        filtered = filtered[filtered["has_stop"] == True]  # noqa: E712

    filtered = filtered.sort_values(["te_id", "aa_len"], ascending=[True, False])
    filtered = filtered.groupby("te_id", as_index=False).head(int(max_orfs_per_te))
    filtered.to_csv(output_table, sep="\t", index=False)

    keep_ids = set(filtered["orf_id"].tolist())
    nt_records = load_fasta_dict(input_nt_fasta)
    aa_records = load_fasta_dict(input_aa_fasta)
    write_subset_fasta(nt_records, keep_ids, output_nt_fasta)
    write_subset_fasta(aa_records, keep_ids, output_aa_fasta)

    logger.info("Retained %d ORFs after filtering", len(filtered))


if __name__ == "__main__":
    try:
        main(
            snakemake.input.table,
            snakemake.input.nt_fasta,
            snakemake.input.aa_fasta,
            snakemake.output.table,
            snakemake.output.nt_fasta,
            snakemake.output.aa_fasta,
            snakemake.params.min_orf_aa,
            snakemake.params.max_orfs_per_te,
            snakemake.params.require_stop_codon,
            snakemake.log[0] if snakemake.log else None,
        )
    except NameError:
        import argparse

        parser = argparse.ArgumentParser(description="Filter ORFs by length and completeness")
        parser.add_argument("--input-table", required=True)
        parser.add_argument("--input-nt-fasta", required=True)
        parser.add_argument("--input-aa-fasta", required=True)
        parser.add_argument("--output-table", required=True)
        parser.add_argument("--output-nt-fasta", required=True)
        parser.add_argument("--output-aa-fasta", required=True)
        parser.add_argument("--min-orf-aa", type=int, default=100)
        parser.add_argument("--max-orfs-per-te", type=int, default=3)
        parser.add_argument("--require-stop-codon", action=argparse.BooleanOptionalAction, default=True)
        parser.add_argument("--log-file")
        args = parser.parse_args()
        main(
            args.input_table,
            args.input_nt_fasta,
            args.input_aa_fasta,
            args.output_table,
            args.output_nt_fasta,
            args.output_aa_fasta,
            args.min_orf_aa,
            args.max_orfs_per_te,
            args.require_stop_codon,
            args.log_file,
        )
