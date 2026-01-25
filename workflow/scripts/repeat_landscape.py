import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import gzip
import os


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

    # Read with flexible space delimiter, handle multi-space issues
    df = pd.read_csv(
        filepath, 
        skiprows=3, 
        sep=r"\s+", 
        names=cols, 
        engine="python",
        on_bad_lines="skip"
    )

    # Calculate length of each hit
    df["length"] = df["q_end"] - df["q_start"]

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

    sns.barplot(
        data=plot_data, x="div_bin", y="genome_perc", hue="class_family", dodge=False
    )

    plt.title("Genomic Repeat Landscape", fontsize=16)
    plt.xlabel("Divergence (Kimura distance %)", fontsize=12)
    plt.ylabel("Percentage of Genome", fontsize=12)
    plt.legend(title="TE Family", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(output_img, dpi=300)
    print(f"Plot saved to {output_img}")


def main(genome: str, out_dataframe: str, landscape: str) -> None:
    """
    Main execution function to generate a repeat landscape.

    Args:
        genome: Path to the genome FASTA file.
        out_dataframe: Path to the RepeatMasker .out file.
        landscape: Path to save the output plot image.
    """
    genome_size = get_genome_size(genome)
    df = parse_repeatmasker_out(out_dataframe, genome_size)
    plot_landscape(df, landscape)


if __name__ == "__main__":
    main(
        snakemake.input.genome, snakemake.repeatmasker_out, snakemake.output.plot
    )
