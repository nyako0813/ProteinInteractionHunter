"""Deterministic flat candidate and warning summary writers."""

import csv
import io
from collections import Counter
from collections.abc import Sequence
from pathlib import Path
from typing import Any

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
)


def _list(values: list[str]) -> str:
    return "|".join(sorted(values))


def _value(value: Any) -> Any:
    return "" if value is None else value


class CandidateTableTsvWriter:
    def write(self, run_id: str, candidates: Sequence[CandidateProtein], path: Path) -> Path:
        output_path = path.expanduser().resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        buffer = io.StringIO(newline="")
        writer = csv.DictWriter(
            buffer, fieldnames=CANDIDATE_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for candidate in sorted(candidates, key=lambda item: (item.query_id, item.protein_id)):
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
