"""Generate or validate checked-in JSON Schemas against Pydantic models."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from jsonschema import Draft202012Validator  # type: ignore[import-untyped]
from pydantic import BaseModel

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from protein_interaction_hunter.models.evidence import CandidateEvidenceBundle  # noqa: E402
from protein_interaction_hunter.models.run import RunManifest  # noqa: E402
from protein_interaction_hunter.models.structure_queue import (  # noqa: E402
    StructurePredictionQueueEntry,
)
from protein_interaction_hunter.schemas.versions import (  # noqa: E402
    SCHEMA_VERSIONS,
    SchemaName,
)

SCHEMA_DIRECTORY = ROOT / "schemas"
SCHEMA_FILES: dict[SchemaName, tuple[type[BaseModel], str]] = {
    SchemaName.CANDIDATE_EVIDENCE_BUNDLE: (
        CandidateEvidenceBundle,
        "candidate_evidence_bundle.schema.json",
    ),
    SchemaName.RUN_MANIFEST: (RunManifest, "run_manifest.schema.json"),
    SchemaName.STRUCTURE_PREDICTION_QUEUE: (
        StructurePredictionQueueEntry,
        "structure_prediction_queue.schema.json",
    ),
}


def generated_schema(name: SchemaName) -> dict[str, Any]:
    model, _ = SCHEMA_FILES[name]
    schema = model.model_json_schema(mode="validation")
    schema["$schema"] = "https://json-schema.org/draft/2020-12/schema"
    schema["$id"] = (
        f"https://schemas.protein-interaction-hunter.local/{name.value}/{SCHEMA_VERSIONS[name]}"
    )
    return schema


def schema_text(name: SchemaName) -> str:
    return json.dumps(generated_schema(name), ensure_ascii=False, indent=2, sort_keys=True) + "\n"


def write_schemas() -> None:
    SCHEMA_DIRECTORY.mkdir(parents=True, exist_ok=True)
    for name, (_, filename) in SCHEMA_FILES.items():
        (SCHEMA_DIRECTORY / filename).write_text(schema_text(name), encoding="utf-8", newline="\n")


def validate_schemas() -> list[str]:
    failures: list[str] = []
    for name, (_, filename) in SCHEMA_FILES.items():
        path = SCHEMA_DIRECTORY / filename
        if not path.is_file():
            failures.append(f"missing checked-in schema: {path}")
            continue
        if path.read_text(encoding="utf-8") != schema_text(name):
            failures.append(f"schema drift detected: {path}")
            continue
        try:
            schema = json.loads(path.read_text(encoding="utf-8"))
            Draft202012Validator.check_schema(schema)
        except Exception as exc:
            failures.append(f"invalid Draft 2020-12 schema {path}: {exc}")
    return failures


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write",
        action="store_true",
        help="Regenerate checked-in schemas from the Pydantic models.",
    )
    args = parser.parse_args()
    if args.write:
        write_schemas()
        print(f"Wrote {len(SCHEMA_FILES)} schemas.")
        return 0
    failures = validate_schemas()
    if failures:
        for failure in failures:
            print(failure, file=sys.stderr)
        return 1
    print(f"Validated {len(SCHEMA_FILES)} schemas with no drift.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
