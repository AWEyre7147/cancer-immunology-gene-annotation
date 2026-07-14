"""Command-line entry point for the gene annotation pipeline."""
import argparse
from .pipeline import run_pipeline

def main():
    # Configure the command-line interface.
    # Each argument corresponds to an option that can be supplied when
    # running the pipeline from a terminal.
    parser = argparse.ArgumentParser(
        description = "Annotate human or mouse gene symbols from a CSV file."
    )

    # Required arguments.
    parser.add_argument(
        "input_csv",
        help = "Path to the input CSV containing gene symbols."
    )

    parser.add_argument(
        "--organism",
        required = True,
        choices = ["human", "mouse"],
        help = "Organism represented by the input gene symbols."
    )

    # Optional arguments.
    parser.add_argument(
        "--output-dir",
        default = "results",
        help = "Directory where annotation tables will be written."
    )

    parser.add_argument(
        "--gene-column",
        help = "Column containing gene symbols. If omitted, the pipeline "
             "attempts to detect the column automatically."
    )

    parser.add_argument(
        "--annotation-dir",
        default = "data/annotation",
        help = "Directory containing local annotation databases."
    )

    parser.add_argument(
        "--skip-remote",
        action = "store_true",
        help = "Skip remote annotation services and use only local resources."
    )

    # Parse the supplied command-line arguments.
    args = parser.parse_args()

    # Execute the annotation pipeline.
    tables = run_pipeline(
        input_csv = args.input_csv,
        organism = args.organism,
        output_dir = args.output_dir,
        gene_column = args.gene_column,
        annotation_dir = args.annotation_dir,
        use_remote = not args.skip_remote,
    )

    # Print a brief summary of the generated output tables.
    print(
        "Wrote:",
        ", ".join(
            f"{name} ({len(table)} rows)"
            for name, table in tables.items()
        ),
    )


# Execute the command-line interface when this module is run directly.
if __name__ == "__main__":
    main()