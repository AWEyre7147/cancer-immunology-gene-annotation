"""Adapters for versioned local annotation tables."""
import pandas as pd

def annotate_hgnc(
    human_genes: list[str],
    hgnc: pd.DataFrame,
) -> pd.DataFrame:
    
    """
    Annotate human genes using the HGNC annotation table.

    Parameters
    ----------
    human_genes : list of str
        Human gene symbols to annotate.

    hgnc : pandas.DataFrame
        HGNC annotation table.

    Returns
    -------
    pandas.DataFrame
        HGNC annotations for the supplied genes.

    Notes
    -----
    Gene symbols are matched in a case-insensitive manner while preserving
    the original order of the input gene list.
    """

    # Store the original input order so the output table can be returned
    # in the same order as the user's input.
    genes = pd.DataFrame({"input_gene": human_genes,
                          "input_order": range(len(human_genes))}
                        )

    # Create case-insensitive lookup keys for the input genes.
    genes["key"] = genes["input_gene"].str.upper()

    # Create matching lookup keys for the HGNC annotation table.
    data = hgnc.copy()
    data["key"] = data["Approved symbol"].astype(str).str.upper()

    # Merge the input genes with the HGNC annotations.
    # A left join ensures that every input gene is retained even if no
    # HGNC annotation is available.
    result = (genes.merge(data, on = "key", how = "left")
        .sort_values("input_order")
        .drop(columns = ["key", "input_order"])
             )

    # Return the annotated gene table.
    return result


def annotate_ligand_receptor(
    human_genes: list[str],
    lr: pd.DataFrame,
) -> pd.DataFrame:
    
    """
    Identify ligand-receptor interactions involving the supplied genes.

    Parameters
    ----------
    human_genes : list of str
        Human gene symbols to search for.

    lr : pandas.DataFrame
        Ligand-receptor interaction table.

    Returns
    -------
    pandas.DataFrame
        Ligand-receptor interactions involving the supplied genes.

    Notes
    -----
    A gene may appear as either the ligand (source) or receptor (target)
    in an interaction. Both roles are searched independently and returned
    in a single table.
    """

    # Build a case-insensitive lookup while preserving the original
    # capitalization of the input gene symbols.
    keys = {gene.upper(): gene
            for gene in human_genes
           }

    # Work on a copy of the annotation table.
    data = lr.copy()

    # Create case-insensitive versions of the ligand and receptor
    # gene-symbol columns.
    source = data["source (gene ID)"].astype(str).str.upper()
    target = data["target (gene ID)"].astype(str).str.upper()

    # Identify interactions where the input gene is the ligand.
    source_hits = data[source.isin(keys)].copy()

    source_hits.insert(0,
                       "input_gene",
                       source[source.isin(keys)].map(keys)
                      )

    source_hits.insert(1,
                       "lr_role",
                       "source"
                      )

    # Identify interactions where the input gene is the receptor.
    target_hits = data[target.isin(keys)].copy()

    target_hits.insert(0,
                       "input_gene",
                       target[target.isin(keys)].map(keys)
                      )

    target_hits.insert(1,
                       "lr_role",
                       "target"
                       )

    # Combine ligand and receptor matches into a single annotation table.
    return pd.concat([source_hits, target_hits],
                     ignore_index = True
                    )