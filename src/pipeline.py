"""Pipeline orchestration: coordinates adapters without hiding their outputs."""

from pathlib import Path
import time
import pandas as pd

from .annotate_local_files import annotate_hgnc, annotate_ligand_receptor
from .input_output import read_gene_csv, write_tables
from .organism_calls import organism_name
from .remote import (
    gprofiler_convert,
    gprofiler_enrichment,
    mouse_to_human,
    mygene_annotations,
)

def _find_first_existing(
    directory: Path,
    filenames: list[str],
) -> Path | None:
    
    """
    Find the first supported annotation file present in a directory.

    Parameters
    ----------
    directory : pathlib.Path
        Directory containing local annotation files.

    filenames : list of str
        Supported filenames listed in order of preference.

    Returns
    -------
    pathlib.Path or None
        Path to the first existing file, or None if no supported filename
        is found.

    Notes
    -----
    Supporting multiple filenames allows the pipeline to remain compatible
    with annotation files that were downloaded or saved under different
    names.

    The order of ``filenames`` determines which file is selected when more
    than one supported file is present.
    """

    # Check each supported filename in order.
    # Return immediately when the first matching file is found.
    for filename in filenames:
        candidate = directory / filename

        if candidate.exists():
            return candidate

    # Return None when none of the supported files are available.
    return None


def _run_remote_step(
    name,
    function,
    tables,
    errors,
    retries = 2,
    retry_delay = 1.0,
):
    
    """
    Run a remote annotation step with retry and error handling.

    Parameters
    ----------
    name : str
        Name used to store the returned table and identify any error.

    function : callable
        Function that performs the remote request and returns a
        pandas DataFrame.

    tables : dict of str to pandas.DataFrame
        Dictionary receiving the result table.

    errors : list of dict
        List receiving non-fatal pipeline errors.

    retries : int, default = 2
        Number of additional attempts made after the initial request fails.

    retry_delay : float, default = 1.0
        Number of seconds to wait between failed attempts.

    Returns
    -------
    None
        Results are added directly to ``tables``. Failures are added directly
        to ``errors``.

    Raises
    ------
    ValueError
        May be raised if ``retries`` is less than zero.

    Notes
    -----
    Remote biological services may fail temporarily because of network
    interruptions, rate limits, maintenance, or server errors. Retrying the
    request can recover from short-lived failures.

    If every attempt fails, an empty DataFrame is stored under ``name`` so
    the expected output key remains available. The final exception message
    is also recorded in the pipeline error table.
    """

    # Validate the retry setting.
    if retries < 0:
        raise ValueError("retries must be zero or greater")

    # Preserve the most recent exception for the final error report.
    last_exception = None

    # Attempt the request once, followed by the requested number of retries.
    for attempt in range(retries + 1):
        try:
            tables[name] = function()
            return

        except Exception as exc:
            # Remote clients may raise different exception types depending
            # on whether the failure comes from HTTP, networking, or parsing.
            last_exception = exc

            # Wait briefly before another attempt, but do not delay after
            # the final failed request.
            if attempt < retries:
                time.sleep(retry_delay)

    # Preserve a consistent output key after all attempts fail.
    tables[name] = pd.DataFrame()

    # Record the final failure without stopping the remaining pipeline.
    errors.append({"step": name,
                   "error": str(last_exception)}
                  )


def run_pipeline(
    input_csv,
    organism,
    output_dir,
    gene_column=None,
    annotation_dir=None,
    use_remote=True,
):
    """
    Run the gene annotation pipeline and save the resulting tables.

    Parameters
    ----------
    input_csv : str or pathlib.Path
        Path to a CSV file containing gene identifiers.

    organism : str
        Organism represented by the input genes.

        Supported options:
            human
            mouse

    output_dir : str or pathlib.Path
        Directory where annotation tables will be saved as CSV files.

    gene_column : str, optional
        Name of the column containing gene identifiers.

        If None, the input-reading function attempts to identify the gene
        column automatically.

    annotation_dir : str or pathlib.Path, optional
        Directory containing local annotation files.

        Supported HGNC filenames:
            gene_annotation_hgnc.txt
            gene annotation.txt
            hgnc_complete_set.txt

        Expected ligand-receptor filename:
            ligand_receptor.csv

        These local resources contain human annotations. For mouse input,
        mouse genes must first be converted to human orthologs before the
        local annotation files can be queried.

    use_remote : bool, default = True
        Whether to query remote annotation services.

        When True, the pipeline attempts to retrieve:
            g:Profiler identifier conversions
            g:Profiler enrichment results
            MyGene.info annotations
            Mouse-to-human orthologs for mouse input

        Remote-service failures are recorded in the ``run_errors`` table
        rather than stopping the entire pipeline.

    Returns
    -------
    dict of str to pandas.DataFrame
        Dictionary containing all annotation tables generated during the run.

        Possible tables include:
            input_genes
            gprofiler_ids
            gprofiler_enrichment
            mygene
            mouse_to_human_orthologs
            hgnc
            ligand_receptor
            run_errors

        Remote and local table keys are preserved with empty DataFrames when
        a requested step cannot be completed.

    Raises
    ------
    ValueError
        If ``organism`` is not "human" or "mouse".

    OSError
        If input files cannot be read, the output directory cannot be
        created, or output tables cannot be written.

    Notes
    -----
    Each returned table is written to a separate CSV file in ``output_dir``.

    Remote annotation steps are run independently. A failure in one remote
    service does not prevent other remote or local annotations from running.

    For human input, the original gene symbols are used when querying HGNC
    and ligand-receptor resources.

    For mouse input, human orthologs returned by g:Profiler are used for
    HGNC and ligand-receptor annotation. If remote services are disabled or
    ortholog conversion fails, those local outputs remain empty and the
    reason is recorded in ``run_errors``.
    """

    # Validate and normalize the requested organism.
    # Lowercasing allows values such as "Human" and "MOUSE" to be accepted
    # while preserving a small, explicit set of supported organisms.
    organism = organism.lower()
    if organism not in {"human", "mouse"}:
        raise ValueError("organism must be 'human' or 'mouse'")

    # Read and preserve the cleaned input gene list.
    # The original input is stored as an output table so every run records
    # exactly which genes entered the pipeline.
    genes = read_gene_csv(input_csv,
                          gene_column)

    tables = {"input_genes": pd.DataFrame({"input_gene": genes})}

    # Collect errors that should not terminate the entire pipeline.
    # Each error is later written to the run_errors.csv output file.
    errors = []

    # Prepare genes for human-only local resources.
    # Human genes can be queried directly. Mouse genes require an ortholog
    # conversion before they can be matched against HGNC or the supplied
    # ligand-receptor database.
    human_genes = genes if organism == "human" else []

    # Query remote annotation services when enabled.
    # Each request is handled independently and retried before being marked
    # as failed.
    if use_remote:
        service_organism = organism_name(
            organism,
            "gprofiler"
        )

        ### Convert input symbols to standardized gene identifiers.
        _run_remote_step(
            "gprofiler_ids",
            lambda: gprofiler_convert(
                genes,
                service_organism
            ),
            tables,
            errors
        )

        ### Perform list-level functional enrichment analysis.
        _run_remote_step(
            "gprofiler_enrichment",
            lambda: gprofiler_enrichment(
                genes,
                service_organism
            ),
            tables,
            errors
        )

        # Retrieve gene-level annotations from MyGene.info.
        _run_remote_step(
            "mygene",
            lambda: mygene_annotations(
                genes,
                organism_name(
                    organism,
                    "mygene"
                ),
            ),
            tables,
            errors
        )

        # Convert mouse genes to human orthologs when needed.
        # The converted symbols are used only for human-specific local
        # resources; the original mouse list remains unchanged.
        if organism == "mouse":
            _run_remote_step(
                "mouse_to_human_orthologs",
                lambda: mouse_to_human(genes),
                tables,
                errors
            )

            orthologs = tables["mouse_to_human_orthologs"]

            # Extract a unique, non-empty list of human ortholog symbols.
            if not orthologs.empty:
                human_genes = (
                    orthologs.get(
                        "human_gene_name",
                        pd.Series(dtype=str),
                    )
                    .dropna()
                    .astype(str)
                    .str.strip()
                    .loc[lambda values: values.ne("")]
                    .drop_duplicates()
                    .tolist()
                )

    # Locate and query local human annotation resources.
    # Output keys are created consistently even if a file is missing or the
    # pipeline has no human symbols available for matching.
    if annotation_dir:
        annotation_path = Path(annotation_dir)

        # Locate the first supported HGNC annotation file.
        hgnc_path = _find_first_existing(
            annotation_path,
            ["gene_annotation_hgnc.txt",
             "gene annotation.txt",
             "hgnc_complete_set.txt"]
        )

        # Locate the ligand-receptor interaction database.
        ligand_receptor_path = _find_first_existing(
            annotation_path,
            ["ligand_receptor.csv"]
        )

        ### Annotate human symbols with the local HGNC table.
        if human_genes and hgnc_path is not None:
            hgnc_data = pd.read_csv(hgnc_path,
                                    sep = "\t",
                                    low_memory = False
                                   )

            tables["hgnc"] = annotate_hgnc(human_genes,
                                           hgnc_data
                                          )

        else:
            # Preserve the output key so users can distinguish a skipped or
            # failed step from a feature that was never implemented.
            tables["hgnc"] = pd.DataFrame()

            if hgnc_path is None:
                errors.append({"step": "hgnc",
                               "error": ("No supported HGNC annotation file was found.")
                               })

            elif organism == "mouse" and not human_genes:
                errors.append({"step": "hgnc",
                               "error": ("Mouse HGNC annotation requires human orthologs. "
                                         "Enable remote services or provide human input "
                                         "genes.")
                              })

        # Identify ligand-receptor interactions involving human symbols.
        if human_genes and ligand_receptor_path is not None:
            ligand_receptor_data = pd.read_csv(ligand_receptor_path,
                                               low_memory = False
                                              )

            tables["ligand_receptor"] = annotate_ligand_receptor(human_genes,
                                                                 ligand_receptor_data
                                                                )

        else:
            # Preserve an empty table when the resource cannot be queried.
            tables["ligand_receptor"] = pd.DataFrame()

            if ligand_receptor_path is None:
                errors.append({"step": "ligand_receptor",
                               "error": "ligand_receptor.csv was not found."
                              })

            elif organism == "mouse" and not human_genes:
                errors.append({"step": "ligand_receptor",
                               "error": ("Mouse ligand-receptor annotation requires human "
                                         "orthologs. Enable remote services or provide "
                                         "human input genes.")
                              })

    # Convert collected error records into a consistent output table.
    # Supplying the columns explicitly preserves the table structure even
    # when no errors occurred.
    tables["run_errors"] = pd.DataFrame(errors,
                                        columns = ["step", "error"]
                                       )

    # Write every generated table to its own CSV file.
    write_tables(tables,
                 output_dir
                )

    # Return the DataFrames for notebooks and other Python workflows.
    return tables
