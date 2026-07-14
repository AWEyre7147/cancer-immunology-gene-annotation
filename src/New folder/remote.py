"""Thin adapters around remote services; easy to replace or extend."""
def gprofiler_convert(
    genes, 
    organism, 
    target_namespace = "ENSG"
):
    
    """
    Convert gene identifiers using g:Profiler.

    Parameters
    ----------
    genes : list-like
        List of gene identifiers (e.g., gene symbols).

    organism : str
        Organism name recognized by g:Profiler.

    target_namespace : str, default = "ENSG"
        Identifier type to convert the input genes to.

        Common options:
            ENSG                Ensembl Gene ID
            ENST                Ensembl Transcript ID
            ENTREZGENE          NCBI Entrez Gene ID
            UNIPROTSWISSPROT    Reviewed UniProt accession

    Returns
    -------
    pandas.DataFrame
        Table containing the converted gene identifiers returned by
        g:Profiler.

    Raises
    ------
    Exception
        Propagates any exceptions raised by the g:Profiler API or network
        connection.

    Notes
    -----
    This function is a lightweight wrapper around the g:Profiler
    ``convert()`` method. It standardizes identifier conversion within the
    annotation pipeline and returns the results as a pandas DataFrame.
    """
    # Import the g:Profiler client only when this function is called.
    # Delaying the import reduces startup time and isolates optional
    # dependencies to the functions that require them.
    from gprofiler import GProfiler

    # Convert the supplied gene identifiers to the requested namespace.
    return GProfiler(return_dataframe = True).convert(
        organism = organism,
        query = genes,
        target_namespace = target_namespace
    )

def mouse_to_human(
    genes
):

    """
    Convert mouse gene symbols to predicted human orthologs using
    g:Profiler.

    Parameters
    ----------
    genes : list-like
        Mouse gene symbols.

    Returns
    -------
    pandas.DataFrame
        Table linking mouse genes to their predicted human orthologs.

        Columns include:
            mouse_gene
            human_gene
            human_gene_name
            human_gene_description
            human_ensembl_gene_id

    Raises
    ------
    Exception
        Propagates any exceptions raised by the g:Profiler API or network
        connection.

    Notes
    -----
    Mouse genes are converted using the g:Profiler orthology service with
    mouse (mmusculus) as the source organism and human (hsapiens) as the
    target organism.

    This function is intended for annotation resources that are available
    only for human genes, such as HGNC and the ligand-receptor database.
    """
    # Import the g:Profiler client.
    from gprofiler import GProfiler

    # Query g:Profiler for mouse-to-human ortholog relationships.
    # These human orthologs allow mouse datasets to be annotated using
    # resources that are available only for human genes.
    df = GProfiler(return_dataframe = True).orth(
        organism = "mmusculus",
        query = genes,
        target = "hsapiens"
    )

    # Rename the returned columns to descriptive names used throughout
    # the annotation pipeline.
    return df.rename(
        columns={
            "input": "mouse_gene",
            "converted": "human_gene",
            "name": "human_gene_name",
            "description": "human_gene_description",
            "ortholog_ensg": "human_ensembl_gene_id"
        }
    )

def mygene_annotations(
    genes, 
    species
):
    
    """
    Retrieve gene annotations from MyGene.info.

    Parameters
    ----------
    genes : list-like
        List of gene symbols.

    species : str
        Species recognized by MyGene.info.

        Common options:
            human
            mouse

    Returns
    -------
    pandas.DataFrame
        Table containing gene annotations returned by MyGene.info.

    Raises
    ------
    Exception
        Propagates any exceptions raised by the MyGene.info service or
        network connection.

    Notes
    -----
    This function queries MyGene.info using gene symbols and returns a
    standardized set of commonly used annotations, including gene names,
    Entrez Gene IDs, Ensembl Gene IDs, aliases, functional summaries,
    UniProt accessions, Gene Ontology terms, and Reactome and KEGG pathway
    annotations.

    Additional annotation fields can be included by modifying the
    ``fields`` argument in the query.
    """
    # Import the MyGene.info client.
    import mygene

    # Create a client for communicating with the MyGene.info API.
    mg = mygene.MyGeneInfo()

    # Retrieve annotations for each gene symbol.
    # Only a curated subset of commonly used annotation fields is requested
    # to keep the output concise while remaining broadly useful.
    return mg.querymany(
        genes,
        scopes = "symbol",
        species = species,
        fields = [
            "symbol",
            "name",
            "entrezgene",
            "ensembl.gene",
            "alias",
            "summary",
            "uniprot.Swiss-Prot",
            "go.BP.term",
            "go.MF.term",
            "go.CC.term",
            "pathway.reactome",
            "pathway.kegg"
        ],
        as_dataframe = True,
        df_index = False,
    )
