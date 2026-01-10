import subprocess
import os
import shutil
import glob
import sys


def run_repeatmodeler():
    # Access snakemake object
    # In snakemake scripts, the 'snakemake' object is available globally.

    database = snakemake.params.db_basename
    threads = snakemake.threads
    output_file = snakemake.output.families
    log_file = snakemake.log[0]

    command = [
        "RepeatModeler",
        "-threads",
        str(threads),
        "-database",
        database,
    ]

    print(f"Running command: {' '.join(command)}")

    with open(log_file, "w") as log:
        process = subprocess.Popen(
            command, stdout=log, stderr=subprocess.STDOUT, universal_newlines=True
        )
        process.wait()

    if process.returncode != 0:
        print(f"RepeatModeler failed with return code {process.returncode}")
        # Note: In a snakemake script, exiting with non-zero will fail the job.
        sys.exit(process.returncode)

    # Find the output directory (RM_...)
    rm_dirs = glob.glob("RM_*")
    if not rm_dirs:
        print("Error: No RM_ directory found after RepeatModeler run.")
        sys.exit(1)

    target_dir = max(rm_dirs, key=os.path.getmtime)
    print(f"Found RepeatModeler output directory: {target_dir}")

    source_file = os.path.join(target_dir, "consensi.fa.classified")
    if not os.path.exists(source_file):
        print(f"Error: Output file {source_file} not found.")
        sys.exit(1)

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Move the file to the desired output path
    shutil.copy2(source_file, output_file)
    print(f"Successfully copied {source_file} to {output_file}")

    # Cleanup
    shutil.rmtree(target_dir)
    print(f"Deleted directory: {target_dir}")


if __name__ == "__main__":
    run_repeatmodeler()
