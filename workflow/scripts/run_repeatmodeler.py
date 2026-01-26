import subprocess
import os
import shutil
import glob
import sys
import logging
from utils import setup_logging

logger = logging.getLogger(__name__)


def run_repeatmodeler_logic(database, threads, output_file, log_file: str):
    setup_logging(log_file, __name__)

    command = [
        "RepeatModeler",
        "-threads",
        str(threads),
        "-database",
        database,
    ]

    logger.info(f"Running command: {' '.join(command)}")

    # Open in append mode because setup_logging might have already opened it
    with open(log_file, "a") as log:
        process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, universal_newlines=True)
        process.wait()

    if process.returncode != 0:
        logger.error(f"RepeatModeler failed with return code {process.returncode}")
        # Note: In a snakemake script, exiting with non-zero will fail the job.
        sys.exit(process.returncode)

    # Find the output directory (RM_...)
    rm_dirs = glob.glob("RM_*")
    if not rm_dirs:
        logger.error("No RM_ directory found after RepeatModeler run.")
        sys.exit(1)

    target_dir = max(rm_dirs, key=os.path.getmtime)
    logger.info(f"Found RepeatModeler output directory: {target_dir}")

    source_file = os.path.join(target_dir, "consensi.fa.classified")
    if not os.path.exists(source_file):
        logger.error(f"Output file {source_file} not found.")
        sys.exit(1)

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Move the file to the desired output path
    shutil.copy2(source_file, output_file)
    logger.info(f"Successfully copied {source_file} to {output_file}")

    # Cleanup
    shutil.rmtree(target_dir)
    logger.info(f"Deleted directory: {target_dir}")


if __name__ == "__main__":
    try:
        run_repeatmodeler_logic(
            database=snakemake.params.db_basename,
            threads=snakemake.threads,
            output_file=snakemake.output.families,
            log_file=snakemake.log[0],
        )
    except NameError:
        import argparse

        parser = argparse.ArgumentParser(description="Run RepeatModeler")
        parser.add_argument("--database", required=True, help="Path to the database file")
        parser.add_argument("--threads", type=int, default=1, help="Number of threads to use")
        parser.add_argument("--output-file", required=True, help="Path to the output file")
        parser.add_argument("--log-file", required=True, help="Path to the log file")
        args = parser.parse_args()
        run_repeatmodeler_logic(
            database=args.database,
            threads=args.threads,
            output_file=args.output_file,
            log_file=args.log_file,
        )
