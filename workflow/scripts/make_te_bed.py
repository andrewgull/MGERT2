import logging
import os

import pandas as pd
from utils import setup_logging

logger = logging.getLogger(__name__)


def parse_repeatmasker_out(filepath: str) -> pd.DataFrame:
    """
    Parses a RepeatMasker .out file and returns a DataFrame.

    The RepeatMasker .out file is expected to be in a fixed-width format.
    This function skips the first three header lines and returns all hits.

    Args:
        filepath: Path to the RepeatMasker .out file.

    Returns:
        A pandas DataFrame containing RepeatMasker results.
    """
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

    def on_bad_line(bad_line):
        logger.warning(f"Skipping bad line in {filepath}: {bad_line}")
        return None

    # Read with flexible space delimiter, handle multi-space issues
    df = pd.read_csv(
        filepath,
        skiprows=3,
        sep=r"\s+",
        names=cols,
        engine="python",
        on_bad_lines=on_bad_line,
    )

    return df


def filter_by_te_name(df: pd.DataFrame, te_name: str) -> pd.DataFrame:
    """
    Filters a RepeatMasker DataFrame by TE name.

    Filters rows where the class_family column contains the specified TE name.

    Args:
        df: DataFrame with RepeatMasker results.
        te_name: The TE name to filter by (will match if substring in class_family).

    Returns:
        Filtered DataFrame.
    """
    filtered = df[
        df["class_family"].str.contains(te_name, na=False, case=False, regex=False)
    ]
    logger.info(
        f"Found {len(filtered)} hits for '{te_name}' out of {len(df)} total hits."
    )
    return filtered


def convert_to_bed(df: pd.DataFrame) -> pd.DataFrame:
    """
    Converts RepeatMasker coordinates to BED format.

    RepeatMasker uses 1-based inclusive coordinates.
    BED format uses 0-based half-open intervals (chrom, start, end).

    Strand conversion: '+' stays '+', 'C' (complement) becomes '-'.

    Args:
        df: Filtered RepeatMasker DataFrame.

    Returns:
        DataFrame with BED format columns:
        chrom, start, end, name (repeat_id), score, strand
    """
    bed_df = pd.DataFrame()
    bed_df["chrom"] = df["query"]
    bed_df["start"] = df["q_start"] - 1  # Convert from 1-based to 0-based
    bed_df["end"] = df["q_end"]  # q_end already inclusive in RM, so no adjustment
    bed_df["name"] = "repeat_" + df["id"].astype(str)
    bed_df["score"] = df["sw_score"]
    bed_df["strand"] = df["strand"].replace({"C": "-", "+": "+"})
    unexpected_strands = bed_df["strand"][~bed_df["strand"].isin(["+", "-"])]
    if not unexpected_strands.empty:
        logger.warning(
            f"Found {len(unexpected_strands)} rows with unexpected strand values"
        )

    return bed_df


def write_bed(bed_df: pd.DataFrame, output_path: str) -> None:
    """
    Writes a BED format file.

    Standard BED format: chrom, start, end, name, score, strand (tab-separated).

    Args:
        bed_df: DataFrame with BED format columns.
        output_path: Path to output BED file.
    """
    bed_df.to_csv(
        output_path,
        sep="\t",
        header=False,
        index=False,
        lineterminator="\n",
    )
    logger.info(f"Wrote {len(bed_df)} regions to {output_path}")


def find_out_file(path, sample=None):
    """
    Finds the .out file in a directory.
    """
    if os.path.isfile(path):
        return path

    if os.path.isdir(path):
        # If sample is provided, look for sample.fasta.out
        if sample:
            potential_file = os.path.join(path, f"{sample}.fasta.out")
            if os.path.exists(potential_file):
                return potential_file

        # Otherwise, look for any .out file
        files = [f for f in os.listdir(path) if f.endswith(".out")]
        if len(files) == 1:
            return os.path.join(path, files[0])
        elif len(files) > 1:
            # Prefer {sample}.fasta.out if multiple exist and sample is known
            if sample:
                matching = [f for f in files if f.startswith(sample)]
                if matching:
                    return os.path.join(path, matching[0])
            logger.warning(
                f"Multiple .out files found in {path}, selecting the first one: {files[0]}"
            )
            return os.path.join(path, files[0])

    raise FileNotFoundError(f"Could not find RepeatMasker .out file at {path}")


def main(
    repeatmasker_input: str,
    bed_output: str,
    te_name: str,
    log_file: str,
    sample: str | None = None,
) -> None:
    """
    Main execution logic.

    Args:
        repeatmasker_input: Path to RepeatMasker .out file or directory.
        bed_output: Path to output BED file.
        te_name: TE name to filter by.
        log_file: Optional path to log file.
        sample: Optional sample name to help find the .out file.
    """
    setup_logging(log_file, __name__)

    # Resolve the .out file path if a directory was provided
    repeatmasker_out = find_out_file(repeatmasker_input, sample)

    logger.info(f"Extracting {te_name} sequences from {repeatmasker_out}")

    # Parse and filter
    df = parse_repeatmasker_out(repeatmasker_out)
    filtered_df = filter_by_te_name(df, te_name)

    if len(filtered_df) == 0:
        logger.warning(f"No hits found for '{te_name}'")

    # Convert to BED and write
    bed_df = convert_to_bed(filtered_df)
    write_bed(bed_df, bed_output)
    logger.info(f"Successfully created BED file with {len(bed_df)} regions")


if __name__ == "__main__":
    try:
        main(
            snakemake.input.repeatmasker_dir,
            snakemake.output.bed,
            snakemake.params.te_name,
            snakemake.log[0] if snakemake.log else None,
            sample=snakemake.wildcards.sample,
        )
    except NameError:
        import argparse

        parser = argparse.ArgumentParser(
            description="Convert RepeatMasker .out hits to BED format, filtered by TE name"
        )
        parser.add_argument(
            "--repeatmasker-out",
            required=True,
            help="Path to RepeatMasker .out file or directory",
        )
        parser.add_argument(
            "--bed-output",
            required=True,
            help="Path to output BED file",
        )
        parser.add_argument(
            "--te-name",
            required=True,
            help="TE name to filter by",
        )
        parser.add_argument(
            "--sample",
            help="Optional sample name to help find the .out file if a directory is provided",
        )
        parser.add_argument(
            "--log-file",
            help="Path to log file",
        )
        args = parser.parse_args()
        main(
            args.repeatmasker_out,
            args.bed_output,
            args.te_name,
            args.log_file,
            sample=args.sample,
        )
