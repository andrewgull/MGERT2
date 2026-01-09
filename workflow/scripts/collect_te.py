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
        if record.id.startswith(te_name):
            collection.append(record)
    return collection


def main(input_fasta, output_fasta, te_name):
    """
    Main execution logic.
    """
    collection = collect_te(input_fasta, te_name)
    SeqIO.write(collection, output_fasta, "fasta")


if __name__ == "__main__":
    main(
        snakemake.input.families,
        snakemake.output.collection,
        snakemake.params.te_name
    )