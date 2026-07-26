"""Validated protein-by-species phylogenetic profile observations."""

from typing import Annotated

from pydantic import StringConstraints

from protein_interaction_hunter.models.base import StrictModel

NonEmptyStr = Annotated[str, StringConstraints(strip_whitespace=True, min_length=1)]


class PhylogeneticProfileObservation(StrictModel):
    protein_id: NonEmptyStr
    species_id: NonEmptyStr
    presence: bool | None
    taxonomic_group: NonEmptyStr | None = None
    source: NonEmptyStr | None = None
    source_record_id: NonEmptyStr | None = None
