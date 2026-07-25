"""Deterministic UTF-8 JSON and JSONL writers."""

import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from pydantic import BaseModel

from protein_interaction_hunter.exceptions import SerializationError
from protein_interaction_hunter.models.evidence import CandidateEvidenceBundle
from protein_interaction_hunter.models.run import RunManifest


def deterministic_json(model: BaseModel) -> str:
    """Serialize a model without NaN/Infinity and with stable key ordering."""
    try:
        return json.dumps(
            model.model_dump(mode="json"),
            allow_nan=False,
            ensure_ascii=False,
            separators=(",", ":"),
            sort_keys=True,
        )
    except (TypeError, ValueError) as exc:
        raise SerializationError(f"Could not serialize {type(model).__name__}: {exc}") from exc


class JsonlEvidenceBundleWriter:
    def write(self, records: Sequence[CandidateEvidenceBundle], path: Path) -> Path:
        output_path = path.expanduser().resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        ordered = sorted(records, key=lambda item: (item.run_id, item.query_id, item.candidate_id))
        text = "".join(f"{deterministic_json(record)}\n" for record in ordered)
        output_path.write_text(text, encoding="utf-8", newline="\n")
        return output_path

    def read(self, path: Path) -> list[CandidateEvidenceBundle]:
        records: list[CandidateEvidenceBundle] = []
        for line_number, line in enumerate(
            path.read_text(encoding="utf-8").splitlines(), start=1
        ):
            if not line.strip():
                continue
            try:
                records.append(CandidateEvidenceBundle.model_validate_json(line))
            except ValueError as exc:
                raise SerializationError(
                    f"Invalid evidence JSONL on line {line_number}: {exc}"
                ) from exc
        return records


class JsonRunManifestWriter:
    def write(self, manifest: RunManifest, path: Path) -> Path:
        output_path = path.expanduser().resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        parsed: Any = json.loads(deterministic_json(manifest))
        output_path.write_text(
            json.dumps(parsed, allow_nan=False, ensure_ascii=False, indent=2, sort_keys=True)
            + "\n",
            encoding="utf-8",
            newline="\n",
        )
        return output_path
