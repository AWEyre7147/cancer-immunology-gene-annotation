"""Central organism-name translation for external services."""
# Mapping between user-facing organism names and the identifiers
# expected by external annotation services. New organisms or services
# can be supported by extending this dictionary without modifying the
# rest of the pipeline.

ORGANISMS = {
    "human": {
        "mygene": "human",
        "gprofiler": "hsapiens",
        "enrichr": "human",
    },
    "mouse": {
        "mygene": "mouse",
        "gprofiler": "mmusculus",
        "enrichr": "mouse",
    },
}


def organism_name(
    organism: str,
    service: str,
) -> str:
    
    """
    Return the organism name required by a specific annotation service.

    Parameters
    ----------
    organism : str
        Common organism name.

        Supported options:
            human
            mouse

    service : str
        Annotation service requiring an organism-specific identifier.

        Supported options depend on the ORGANISMS mapping. Examples include:
            gprofiler
            mygene

    Returns
    -------
    str
        Organism identifier formatted for the requested service.

        Examples:
            human + gprofiler -> hsapiens
            mouse + gprofiler -> mmusculus
            human + mygene -> human
            mouse + mygene -> mouse

    Raises
    ------
    ValueError
        If either the organism or service is not present in the
        ORGANISMS mapping.

    Notes
    -----
    This function provides a single location for translating user-friendly
    organism names into the naming conventions expected by external
    annotation services. Centralizing these mappings makes it easy to add
    support for additional organisms or databases without modifying the
    rest of the pipeline.
    """

    # Normalize the inputs so organism and service names are
    # matched in a case-insensitive manner.
    organism = organism.lower()
    service = service.lower()

    # Look up the service-specific organism identifier.
    # Raise a descriptive error if either the organism or service is
    # not supported by the pipeline.
    try:
        return ORGANISMS[organism][service]
    except KeyError as exc:
        raise ValueError(
            f"Unsupported organism/service: {organism}/{service}"
        ) from exc