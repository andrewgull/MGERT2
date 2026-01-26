from Bio import SeqIO


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


def main(input_fasta, output_fasta, te_name):
    """
    Main execution logic.
    """
    collection = collect_te(input_fasta, te_name)
    SeqIO.write(collection, output_fasta, "fasta")


if __name__ == "__main__":
    try:
        main(snakemake.input.families, snakemake.output.collection, snakemake.params.te_name)
    except NameError:
        import argparse
        parser = argparse.ArgumentParser(description="Collect TE sequences from a FASTA file")
        parser.add_argument("--input-fasta", required=True, help="Path to the input FASTA file")
        parser.add_argument("--output-fasta", required=True, help="Path to the output FASTA file")
        parser.add_argument("--te-name", required=True, help="Prefix of the sequence IDs to collect")
        args = parser.parse_args()
        main(args.input_fasta, args.output_fasta, args.te_name)
