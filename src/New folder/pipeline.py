"""Pipeline orchestration: coordinates adapters without hiding their outputs."""
from pathlib import Path
import pandas as pd
from .input_output import read_gene_csv, write_tables
from .organism_calls import organism_name
from .annotate_local_files import annotate_hgnc, annotate_ligand_receptor
from .remote import gprofiler_convert, mouse_to_human, mygene_annotations

def run_pipeline(
    input_csv,
    organism, 
    output_dir, 
    gene_column = None, 
    annotation_dir = None, 
    use_remote = True
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
        Column containing the gene identifiers.
    
        Name of the column containing gene identifiers. If None, the
        input-reading function attempts to determine the
        appropriate column automatically.
    
    annotation_dir : str or pathlib.Path, optional
        Directory containing local annotation files.
    
        Expected files:
            gene_annotation_hgnc.txt
            ligand_receptor.csv
    
        These local resources contain human annotations. For mouse input,
        mouse genes are first converted to human orthologs before the local
        annotation files are queried.
    
    use_remote : bool, default = True
        Whether to query remote annotation services.
    
        When True, the pipeline attempts to retrieve:
            g:Profiler identifier conversions
            MyGene.info annotations
            Mouse-to-human orthologs for mouse input
    
        Remote-service failures are recorded in the run_errors table rather
        than stopping the entire pipeline.
    
    Returns
    -------
    dict of str to pandas.DataFrame
        Dictionary containing the annotation tables generated during the run.
    
        Possible tables include:
            input_genes
            gprofiler_ids
            mygene
            mouse_to_human_orthologs
            hgnc
            ligand_receptor
            run_errors
    
        The available tables depend on the organism, enabled services,
        available annotation files, and whether any annotation steps fail.
    
    Raises
    ------
    ValueError
        If organism is not "human" or "mouse".
    
    Notes
    -----
    Each table is written to a separate CSV file in output_dir.
    
    HGNC and ligand-receptor annotations are only generated when:
        annotation_dir is provided,
        the expected local file exists, and
        at least one human gene symbol is available.
    
    For human input, the original gene list is used for human annotation
    resources. For mouse input, human orthologs returned by g:Profiler are
    used instead.
    """
    # Validate the requested organism.
    # The pipeline currently supports only human and mouse.
    organism = organism.lower()
    if organism not in {"human", "mouse"}:
        raise ValueError("organism must be 'human' or 'mouse'")
    
    # Read the input gene list.
    genes = read_gene_csv(
        input_csv, 
        gene_column
    )
    tables = {
        "input_genes": pd.DataFrame({"input_gene": genes})
    }
    
    # Collect non-fatal errors encountered during annotation.
    errors = []
    
    # Human annotation resources can be queried directly for human input.
    human_genes = genes if organism == "human" else []
    
    # Query remote annotation services.
    # Each service is executed independently so a failure in one service does
    # not prevent the remaining annotations from being generated.
    if use_remote:
        for name, fn in {
            "gprofiler_ids": lambda: gprofiler_convert(
                genes,
                organism_name(organism, "gprofiler"),
            ),
            "mygene": lambda: mygene_annotations(
                genes,
                organism_name(organism, "mygene"),
            ),
        }.items():
            try:
                tables[name] = fn()
            except Exception as exc:
                errors.append(
                    {"step": name, "error": str(exc)}
                )
    
        # Convert mouse genes to predicted human orthologs so they can be
        # matched against annotation resources that contain only human genes.
        if organism == "mouse":
            try:
                orth = mouse_to_human(genes)
                tables["mouse_to_human_orthologs"] = orth
                human_genes = (
                    orth.get("human_gene_name", pd.Series(dtype=str))
                        .dropna()
                        .astype(str)
                        .drop_duplicates()
                        .tolist()
                )
            except Exception as exc:
                errors.append({
                        "step": "mouse_to_human_orthologs",
                        "error": str(exc),
                    }
                )
    
    # Query local annotation resources.
    if annotation_dir and human_genes:
        ann = Path(annotation_dir)
        hgnc_path = ann / "gene_annotation_hgnc.txt"
        lr_path = ann / "ligand_receptor.csv"
    
        # HGNC gene annotation.
        if hgnc_path.exists():
            tables["hgnc"] = annotate_hgnc(
                human_genes,
                pd.read_csv(
                    hgnc_path,
                    sep = "\t",
                    low_memory = False
                )
            )
    
        # Ligand-receptor annotation.
        if lr_path.exists():
            tables["ligand_receptor"] = annotate_ligand_receptor(
                human_genes,
                pd.read_csv(
                    lr_path,
                    low_memory = False
                )
            )
    
    
    # Save any annotation errors encountered during the run.
    tables["run_errors"] = pd.DataFrame(
        errors,
        columns = ["step", "error"]
    )
    
    
    # Write every annotation table to the output directory.
    write_tables(
        tables, 
        output_dir)
    
    # Return all generated tables for interactive use or downstream analysis.
    return tables
