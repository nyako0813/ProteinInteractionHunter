"""Local-only input adapters."""

from protein_interaction_hunter.adapters.local.annotation import LocalAnnotationTsvLoader
from protein_interaction_hunter.adapters.local.fasta import LocalFastaLoader
from protein_interaction_hunter.adapters.local.gff import LocalGff3Loader

__all__ = ["LocalAnnotationTsvLoader", "LocalFastaLoader", "LocalGff3Loader"]
