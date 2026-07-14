"""Input/output helpers kept separate from biological database logic."""
from pathlib import Path
import pandas as pd

def read_gene_csv(
    path: str | Path,
    gene_column: str | None = None,
) -> list[str]:
    
    """
    Read gene symbols from a CSV file.

    Parameters
    ----------
    path : str or pathlib.Path
        Path to the input CSV file.

    gene_column : str, optional
        Name of the column containing gene symbols.

        If None, the function attempts to identify the gene column
        automatically. Automatic detection succeeds only when the CSV
        contains exactly one non-index column.

    Returns
    -------
    list of str
        Unique, non-empty gene symbols in their original order.

    Raises
    ------
    ValueError
        If the gene column cannot be identified automatically.

    ValueError
        If the requested gene column is not present in the CSV file.

    Notes
    -----
    Missing values, empty strings, and duplicate gene symbols are removed.

    Gene symbols are converted to strings and stripped of leading and
    trailing whitespace. Their capitalization is not modified because
    gene-symbol conventions differ between organisms.
    """

    # Read the input CSV into a DataFrame.
    df = pd.read_csv(path)

    # Identify the gene-symbol column when one was not supplied.
    # Columns beginning with "Unnamed" are commonly produced when a
    # DataFrame index is accidentally written to a CSV file. These columns
    # are excluded from automatic gene-column detection.
    if gene_column is None:
        candidates = [
            column
            for column in df.columns
            if not str(column).lower().startswith("unnamed")
        ]

        # Automatic detection is intentionally conservative. Requiring
        # exactly one candidate prevents the pipeline from silently choosing
        # the wrong column in a multi-column input file.
        if len(candidates) != 1:
            raise ValueError(
                "Specify --gene-column when the CSV has multiple data columns."
            )

        gene_column = candidates[0]

    # Confirm that the requested column exists.
    if gene_column not in df.columns:
        raise ValueError(
            f"Gene column {gene_column!r} not found. "
            f"Available columns: {list(df.columns)}"
        )

    # Clean the gene list while preserving its original order.
    # Processing steps:
    #   1. Remove missing values.
    #   2. Convert values to strings.
    #   3. Remove surrounding whitespace.
    #   4. Remove empty strings.
    #   5. Remove duplicate symbols.
    genes = (
        df[gene_column]
        .dropna()
        .astype(str)
        .str.strip()
        .loc[lambda values: values.ne("")]
        .drop_duplicates()
        .tolist()
    )

    # Return the cleaned gene-symbol list.
    return genes


def write_tables(
    tables: dict[str, pd.DataFrame],
    output_dir: str | Path,
) -> list[Path]:
    
    """
    Write annotation tables to separate CSV files.

    Parameters
    ----------
    tables : dict of str to pandas.DataFrame
        Mapping of output table names to DataFrames.

        Each dictionary key becomes the CSV filename. For example,
        a key named ``mygene`` is written as ``mygene.csv``.

    output_dir : str or pathlib.Path
        Directory where the CSV files will be written.

    Returns
    -------
    list of pathlib.Path
        Paths to the CSV files created by the function.

    Raises
    ------
    OSError
        If the output directory cannot be created or a file cannot be
        written.

    AttributeError
        If a value in ``tables`` does not provide a ``to_csv`` method.

    Notes
    -----
    The output directory and any missing parent directories are created
    automatically.

    Existing files with matching names are overwritten.

    DataFrame indexes are not written to the CSV files because they are
    internal pandas row labels rather than biological annotation fields.
    """

    # Create the output directory if it does not already exist.
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Keep a record of every file created during the pipeline run.
    paths: list[Path] = []

    # Write each annotation result as an independent CSV file.
    for name, table in tables.items():
        path = out / f"{name}.csv"

        # Exclude the DataFrame index so only explicit annotation columns
        # appear in the output file.
        table.to_csv(path, index=False)

        paths.append(path)

    # Return the generated file paths for logging or downstream use.
    return paths