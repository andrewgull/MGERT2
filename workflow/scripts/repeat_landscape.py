import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import gzip
import os
import logging
from utils import setup_logging
logger = logging.getLogger(__name__)


def parse_repeatmasker_out(filepath: str, genome_size: int) -> pd.DataFrame:
    """
    Parses a RepeatMasker .out file and calculates the repeat landscape.

    The RepeatMasker .out file is expected to be in a fixed-width format.
    This function skips the first three header lines, calculates the length
    of each hit, groups by divergence (rounded to the nearest integer) and
    repeat class/family, and calculates the percentage of the genome
    occupied by each category.

    Args:
        filepath: Path to the RepeatMasker .out file.
        genome_size: Total size of the genome in base pairs.

    Returns:
        A pandas DataFrame containing 'div_bin', 'class_family', 'length',
        and 'genome_perc'.
    """
    # RepeatMasker .out is a fixed-width format.
    # We skip the first 3 header lines.
    cols = [
        "sw_score",
        "perc_div",
        "perc_del",
        "perc_ins",
        "query",
        "q_start",
        "q_end",
        "q_left",
        "strand",
        "repeat",
        "class_family",
        "r_start",
        "r_end",
        "r_left",
        "id",
    ]
    if genome_size <= 0:
        logger.error(f"Invalid genome_size: {genome_size}. Must be positive.")
        raise ValueError(f"Invalid genome_size: {genome_size}. Must be positive.")

    def on_bad_line(bad_line):
        logger.warning(f"Skipping bad line in {filepath}: {bad_line}")
        return None

    # Read with flexible space delimiter, handle multi-space issues
    df = pd.read_csv(filepath, skiprows=3, sep=r"\s+", names=cols, engine="python", on_bad_lines=on_bad_line)

    # Calculate length of each hit
    df["length"] = df["q_end"] - df["q_start"] + 1

    # We group by divergence (rounded to nearest integer) and Repeat Class
    # Standard classes: LINE, SINE, LTR, DNA, etc.
    df["div_bin"] = df["perc_div"].round()

    # Calculate % of genome for each bin
    landscape = df.groupby(["div_bin", "class_family"])["length"].sum().reset_index()
    landscape["genome_perc"] = (landscape["length"] / genome_size) * 100

    return landscape


def get_genome_size(fasta_path: str) -> int:
    """
    Calculates the total genome size from a FASTA file.

    Sums the number of characters in all sequences, excluding headers
    (lines starting with '>') and newlines. Supports gzipped files.

    Args:
        fasta_path: Path to the FASTA file (can be .gz).

    Returns:
        The total number of bases in the genome.
    """
    size = 0
    open_func = gzip.open if fasta_path.endswith(".gz") else open
    mode = "rt" if fasta_path.endswith(".gz") else "r"

    with open_func(fasta_path, mode) as f:
        for line in f:
            if not line.startswith(">"):
                size += len(line.strip())
    return size


def plot_landscape(data: pd.DataFrame, output_img: str) -> None:
    """
    Generates and saves a Repeat Landscape plot.

    Filters for the top 5 most abundant repeat classes and creates a
    stacked-like bar plot (though using sns.barplot with dodge=False)
    showing genome percentage vs. divergence.

    Args:
        data: DataFrame containing 'div_bin', 'class_family', and 'genome_perc'.
        output_img: Path where the resulting plot will be saved.
    """
    plt.figure(figsize=(12, 7))
    sns.set_style("whitegrid")

    # Filter for major classes to avoid a messy legend
    major_classes = data.groupby("class_family")["genome_perc"].sum().nlargest(5).index
    plot_data = data[data["class_family"].isin(major_classes)]

    sns.barplot(data=plot_data, x="div_bin", y="genome_perc", hue="class_family", dodge=False)

    plt.title("Genomic Repeat Landscape", fontsize=16)
    plt.xlabel("Divergence (Kimura distance %)", fontsize=12)
    plt.ylabel("Percentage of Genome", fontsize=12)
    plt.legend(title="TE Family", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(output_img, dpi=300)
    plt.close()
    logger.info(f"Plot saved to {output_img}")


def main(genome: str, repeatmasker_dir: str, landscape: str, log_file: str = None) -> None:
    """
    Main execution function to generate a repeat landscape.

    Args:
        genome: Path to the genome FASTA file.
        repeatmasker_dir: Directory containing RepeatMasker output files.
        landscape: Path to save the output plot image.
        log_file: Path to the log file.
    """
    setup_logging(log_file, __name__)
    # Find the .out file in the repeatmasker_dir
    out_files = sorted(f for f in os.listdir(repeatmasker_dir) if f.endswith(".out"))
    if not out_files:
        raise FileNotFoundError(f"No .out file found in {repeatmasker_dir}")
    if len(out_files) > 1:
        # Optionally pick one or warn, here we pick the first one
        logger.warning(f"Multiple .out files found in {repeatmasker_dir}. Using {out_files[0]}")

    out_dataframe = os.path.join(repeatmasker_dir, out_files[0])

    genome_size = get_genome_size(genome)
    df = parse_repeatmasker_out(out_dataframe, genome_size)
    plot_landscape(df, landscape)


if __name__ == "__main__":
    try:
        # When run via Snakemake
        main(snakemake.input.genome, snakemake.input.repeatmasker_dir, snakemake.output.plot, snakemake.log[0] if snakemake.log else None)
    except NameError:
        # When run standalone
        import argparse

        parser = argparse.ArgumentParser(description="Generate repeat landscape plot")
        parser.add_argument("--genome", required=True, help="Path to genome FASTA")
        parser.add_argument("--repeatmasker-dir", required=True, help="RepeatMasker output directory")
        parser.add_argument("--output", required=True, help="Output plot path")
        parser.add_argument("--log-file", help="Path to the log file")
        args = parser.parse_args()
        main(args.genome, args.repeatmasker_dir, args.output, args.log_file)
