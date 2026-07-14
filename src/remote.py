"""Thin adapters around remote services; easy to replace or extend."""


def gprofiler_convert(
    genes,
    organism,
    target_namespace="ENSG",
):
    """Convert gene identifiers using g:Profiler."""
    from gprofiler import GProfiler

    return GProfiler(return_dataframe=True).convert(
        organism=organism,
        query=genes,
        target_namespace=target_namespace,
    )


def gprofiler_enrichment(
    genes,
    organism,
):
    """Run functional enrichment analysis for a gene list using g:Profiler."""
    from gprofiler import GProfiler

    return GProfiler(return_dataframe=True).profile(
        organism=organism,
        query=genes,
        user_threshold=0.05,
        significance_threshold_method="g_SCS",
        no_evidences=False,
    )


def mouse_to_human(genes):
    """Convert mouse gene symbols to predicted human orthologs using g:Profiler."""
    from gprofiler import GProfiler

    df = GProfiler(return_dataframe=True).orth(
        organism="mmusculus",
        query=genes,
        target="hsapiens",
    )

    return df.rename(
        columns={
            "input": "mouse_gene",
            "incoming": "mouse_gene",
            "converted": "human_gene",
            "name": "human_gene_name",
            "description": "human_gene_description",
            "ortholog_ensg": "human_ensembl_gene_id",
        }
    )


def mygene_annotations(
    genes,
    species,
):
    """Retrieve gene annotations from MyGene.info."""
    import mygene

    mg = mygene.MyGeneInfo()

    return mg.querymany(
        genes,
        scopes="symbol",
        species=species,
        fields=[
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
            "pathway.kegg",
        ],
        as_dataframe=True,
        df_index=False,
    )
