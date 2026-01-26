import logging
from Bio import SeqIO
from utils import setup_logging

logger = logging.getLogger(__name__)


def collect_te(fasta_file, te_name):
    """
    Collects sequences from a FASTA file that start with a specific TE name.

    Args:
        fasta_file: Path to the input FASTA file or a file-like object.
        te_name: The prefix of the sequence IDs to collect.

    Returns:
        A list of SeqRecord objects.
    """
    records = SeqIO.parse(fasta_file, "fasta")
    collection = []
    for record in records:
        if te_name in record.id:
            collection.append(record)
    return collection


def main(input_fasta, output_fasta, te_name, log_file=None):
    """
    Main execution logic.
    """
    setup_logging(log_file, __name__)
    logger.info(f"Collecting sequences for TE: {te_name}")
    collection = collect_te(input_fasta, te_name)
    logger.info(f"Found {len(collection)} sequences.")
    SeqIO.write(collection, output_fasta, "fasta")


if __name__ == "__main__":
    try:
        main(
            snakemake.input.families,
            snakemake.output.collection,
            snakemake.params.te_name,
            snakemake.log[0] if snakemake.log else None,
        )
    except NameError:
        import argparse

        parser = argparse.ArgumentParser(description="Collect TE sequences from a FASTA file")
        parser.add_argument("--input-fasta", required=True, help="Path to the input FASTA file")
        parser.add_argument("--output-fasta", required=True, help="Path to the output FASTA file")
        parser.add_argument("--te-name", required=True, help="TE name prefix to filter sequences")
        parser.add_argument("--log-file", help="Path to the log file")
        args = parser.parse_args()
        main(args.input_fasta, args.output_fasta, args.te_name, args.log_file)
