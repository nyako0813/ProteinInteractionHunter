"""Deterministic flat candidate and warning summary writers."""

import csv
import io
from collections import Counter
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from protein_interaction_hunter.models.evidence import GenomeContextEvidence
from protein_interaction_hunter.models.protein import CandidateProtein

CANDIDATE_COLUMNS = (
    "run_id",
    "query_id",
    "candidate_id",
    "candidate_description",
    "candidate_disposition",
    "disposition_reasons",
    "sequence_length",
    "gene_id",
    "locus_tag",
    "old_locus_tag",
    "contig",
    "strand",
    "has_coordinate",
    "has_annotation",
    "same_contig_as_query",
    "is_duplicate_sequence",
    "duplicate_group_id",
    "is_fragment_candidate",
    "fragment_reasons",
    "is_hypothetical",
    "identifier_match_status",
    "warnings",
    "same_contig",
    "query_start",
    "query_end",
    "query_strand",
    "candidate_start",
    "candidate_end",
    "candidate_strand",
    "strand_relationship",
    "relative_position",
    "coordinate_position",
    "distance_bp",
    "overlap_bp",
    "intervening_gene_count",
    "intervening_feature_count",
    "feature_index_delta",
    "within_neighborhood_window",
    "context_completeness",
    "gene_context_status",
)


def _list(values: list[str]) -> str:
    return "|".join(sorted(values))


def _value(value: Any) -> Any:
    if value is None:
        return ""
    return value.value if hasattr(value, "value") else value


class CandidateTableTsvWriter:
    def write(
        self,
        run_id: str,
        candidates: Sequence[CandidateProtein],
        path: Path,
        contexts: Mapping[tuple[str, str], GenomeContextEvidence] | None = None,
    ) -> Path:
        output_path = path.expanduser().resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        buffer = io.StringIO(newline="")
        writer = csv.DictWriter(
            buffer, fieldnames=CANDIDATE_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        context_map = contexts or {}
        for candidate in sorted(candidates, key=lambda item: (item.query_id, item.protein_id)):
            context = context_map.get((candidate.query_id, candidate.protein_id))
            writer.writerow(
                {
                    "run_id": run_id,
                    "query_id": candidate.query_id,
                    "candidate_id": candidate.protein_id,
                    "candidate_description": candidate.description,
                    "candidate_disposition": candidate.disposition.value,
                    "disposition_reasons": _list(candidate.disposition_reasons),
                    "sequence_length": candidate.sequence_length,
                    "gene_id": _value(candidate.gene_id),
                    "locus_tag": _value(candidate.locus_tag),
                    "old_locus_tag": _value(candidate.old_locus_tag),
                    "contig": _value(candidate.contig),
                    "strand": _value(candidate.strand),
                    "has_coordinate": candidate.has_coordinate,
                    "has_annotation": candidate.has_annotation,
                    "same_contig_as_query": _value(candidate.same_contig_as_query),
                    "is_duplicate_sequence": candidate.is_duplicate_sequence,
                    "duplicate_group_id": _value(candidate.duplicate_sequence_group),
                    "is_fragment_candidate": candidate.is_fragment_candidate,
                    "fragment_reasons": _list(candidate.fragment_reasons),
                    "is_hypothetical": candidate.is_hypothetical,
                    "identifier_match_status": candidate.identifier_match_status.value,
                    "warnings": _list(candidate.warnings),
                    "same_contig": _value(context.same_contig if context else None),
                    "query_start": _value(context.query_start if context else None),
                    "query_end": _value(context.query_end if context else None),
                    "query_strand": _value(context.query_strand if context else None),
                    "candidate_start": _value(context.candidate_start if context else None),
                    "candidate_end": _value(context.candidate_end if context else None),
                    "candidate_strand": _value(context.candidate_strand if context else None),
                    "strand_relationship": _value(context.strand_relationship if context else None),
                    "relative_position": _value(context.relative_position if context else None),
                    "coordinate_position": _value(context.coordinate_position if context else None),
                    "distance_bp": _value(context.distance_bp if context else None),
                    "overlap_bp": _value(context.overlap_bp if context else None),
                    "intervening_gene_count": _value(
                        context.intervening_gene_count if context else None
                    ),
                    "intervening_feature_count": _value(
                        context.intervening_feature_count if context else None
                    ),
                    "feature_index_delta": _value(context.feature_index_delta if context else None),
                    "within_neighborhood_window": _value(
                        context.within_neighborhood_window if context else None
                    ),
                    "context_completeness": _value(
                        context.context_completeness if context else None
                    ),
                    "gene_context_status": _value(context.status if context else None),
                }
            )
        output_path.write_text(buffer.getvalue(), encoding="utf-8", newline="\n")
        return output_path


class WarningSummaryTsvWriter:
    def write(self, warnings: Sequence[str], path: Path) -> Path:
        output_path = path.expanduser().resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        buffer = io.StringIO(newline="")
        writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
        writer.writerow(("warning", "count"))
        writer.writerows(sorted(Counter(warnings).items()))
        output_path.write_text(buffer.getvalue(), encoding="utf-8", newline="\n")
        return output_path
