"""Public package interface for the gene annotation pipeline."""

# Expose the primary pipeline function at the package level.
# This allows users to import the pipeline directly:
#     from gene_annotation_pipeline import run_pipeline
from .pipeline import run_pipeline

# Define the public objects exported when using:
#     from gene_annotation_pipeline import *
# Only the primary pipeline function is exposed as part of the public API.
__all__ = ["run_pipeline"]

# Current package version.
__version__ = "0.1.0"