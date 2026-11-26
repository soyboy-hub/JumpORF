#!/usr/bin/env python3

__author__ = "Tiurin K."
__author_email__ = "tiurin.kn@gmail.com"

import argparse
import shutil
import subprocess
import sys
import time
from pathlib import Path


MODULE_DIR = Path("/module")


def positive_int(value: str) -> int:
    """Validate that an argument is a positive integer."""
    try:
        value_int = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"expected an integer, got '{value}'"
        )

    if value_int < 1:
        raise argparse.ArgumentTypeError("value must be greater than 0")

    return value_int


def percent_float(value: str) -> float:
    """Validate that an argument is a percentage between 0 and 100."""
    try:
        value_float = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"expected a number, got '{value}'"
        )

    if not 0 <= value_float <= 100:
        raise argparse.ArgumentTypeError(
            "value must be between 0 and 100"
        )

    return value_float


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="repeat_annotation_pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Run the repeat-annotation pipeline: LTR-RTE, DNA-TE, Helitron, "
            "LINE, and aORF detection. Consensus generation is optional."
        ),
        epilog=(''),
    )

    mandatory = parser.add_argument_group("mandatory arguments")

    mandatory.add_argument(
        "--genome",
        type=Path,
        required=True,
        help=(
            "Path to the input genome assembly in FASTA format "
            "(.fasta, .fna, or .fa)."
        ),
    )

    mandatory.add_argument(
        "--outpath",
        type=Path,
        required=True,
        help="Path to the output directory.",
    )

    mandatory.add_argument(
        "--prefix",
        required=True,
        help=(
            "Prefix used for output files."
        ),
    )

    mandatory.add_argument(
        "--processors",
        type=positive_int,
        required=True,
        help="Number of CPU processors/threads to use.",
    )

    mandatory.add_argument(
        "--species",
        required=True,
        help="Species name supported by the installed AUGUSTUS configuration.",
    )

    mandatory.add_argument(
        "--min-identity-to-cluster",
        type=percent_float,
        required=True,
        help=(
            "Minimum sequence identity required to cluster sequences, f.e. 0.8."
        ),
    )

    mandatory.add_argument(
        "--min-size-cluster",
        type=positive_int,
        required=True,
        help=(
            "Minimum number of sequences required to form a cluster."
        ),
    )

    mandatory.add_argument(
        "--min-coverage",
        type=percent_float,
        required=True,
        help=(
            "Minimum alignment coverage required for clustering, f.e. 0.5."
        ),
    )

    mandatory.add_argument(
        "--fallback-mode",
        choices=("yes", "no"),
        required=True,
        help="Enable the DANTE_LTR fallback mode: 'yes' or 'no'.",
    )

    optional = parser.add_argument_group("optional arguments")

    optional.add_argument(
        "--run-consensus",
        action="store_true",
        help=(
            "Run the optional consensus-generation module after detection."
        ),
    )

    optional.add_argument(
        "--overwrite",
        action="store_true",
        help=(
            "Delete and recreate the output directory if it already exists."
        ),
    )

    optional.add_argument(
        "--debug",
        action="store_true",
        help=(
            "Keep the temporary directory after the pipeline finishes."
        ),
    )

    return parser.parse_args()


def run_module(
    script_name: str,
    arguments: list,
) -> None:
    """Run one pipeline module and stop if it fails."""
    script = MODULE_DIR / script_name

    if not script.is_file():
        raise FileNotFoundError(
            f"Required module not found: {script}"
        )

    command = [
        sys.executable,
        str(script),
        *map(str, arguments),
    ]

    print(
        f"INFO: Running: {' '.join(command)}",
        flush=True,
    )

    result = subprocess.run(command)

    if result.returncode != 0:
        raise RuntimeError(
            f"Module failed with exit code {result.returncode}: "
            f"{script_name}"
        )


def prepare_output_directory(
    args: argparse.Namespace,
) -> tuple[Path, Path]:
    """Create or overwrite the output and temporary directories."""
    outdir = args.outpath.resolve()
    outdir_tmp = outdir / "tmp"

    if outdir.exists():
        if not args.overwrite:
            print(
                f"ERROR: Output directory already exists: {outdir}\n"
                "Use --overwrite or choose a new output directory.",
                file=sys.stderr,
            )
            sys.exit(1)

        print(
            f"INFO: Overwriting existing directory: {outdir}",
            flush=True,
        )
        shutil.rmtree(outdir)

    outdir.mkdir(parents=True, exist_ok=True)
    outdir_tmp.mkdir(parents=True, exist_ok=True)

    return outdir, outdir_tmp


def print_banner() -> None:
    print(r"      _                        ___   _____  _____ ", flush=True)
    print(r"     | |                      / _ \ |  __ ||  ___)", flush=True)
    print(
        r"     | |_   _ ________   ___ | | | || |__) | |__"
        r"                       _________",
        flush=True,
    )
    print(
        r" _   | | | | |  _   _ \ / _ \| | | ||  _  /|  __)"
        r"          _______    /     ____\____",
        flush=True,
    )
    print(
        r"| |__| | | | | | | | | | (_) | |_| || | \ \| |"
        r"     _____  /       \  /     |_________|->  ",
        flush=True,
    )
    print(
        r" \____/ \___/|_| |_| |_|  __/ \___/ |_|  \_\_|"
        r"    /     \/         \/",
        flush=True,
    )
    print(r"                       | |", flush=True)
    print(r"                       |_|", flush=True)


def main() -> None:
    args = get_args()

    if not args.genome.is_file():
        print(
            f"ERROR: Genome file does not exist: {args.genome}",
            file=sys.stderr,
        )
        sys.exit(1)

    outdir, outdir_tmp = prepare_output_directory(args)

    print_banner()

    start = time.time()

    common_arguments = [
        args.genome,
        args.processors,
        outdir,
        args.prefix,
        outdir_tmp,
    ]

    try:
        # Run LTR-RTE detection module.
        run_module(
            "run_LTR_RTE_detection.py",
            [
                *common_arguments,
                args.fallback_mode,
            ],
        )

        # Run DNA-TE detection module.
        run_module(
            "run_DNA_TE_detection.py",
            common_arguments,
        )

        # Run Helitron detection module.
        run_module(
            "run_Helitron_detection.py",
            common_arguments,
        )

        # Run LINE detection module.
        run_module(
            "run_LINE_detection.py",
            common_arguments,
        )

        # Run aORF detection module.
        run_module(
            "run_aORF_detection.py",
            [
                *common_arguments,
                args.species,
                args.min_identity_to_cluster,
                args.min_size_cluster,
                args.min_coverage,
            ],
        )

        # Run optional consensus-generation module.
        if args.run_consensus:
            run_module(
                "run_consensus.py",
                common_arguments,
            )
        else:
            print(
                "INFO: Consensus generation skipped. "
                "Use --run-consensus to enable it.",
                flush=True,
            )

    except (FileNotFoundError, RuntimeError) as error:
        print(
            f"ERROR: {error}",
            file=sys.stderr,
        )
        sys.exit(1)

    elapsed = time.time() - start

    print(
        f"INFO: Pipeline finished successfully. "
        f"Total execution time: {elapsed:.2f} seconds.",
        flush=True,
    )

    if not args.debug:
        shutil.rmtree(outdir_tmp)


if __name__ == "__main__":
    main()
    

