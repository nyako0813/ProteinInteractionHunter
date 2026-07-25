# Schemas

All checked-in JSON Schemas use JSON Schema Draft 2020-12 and are generated from the Pydantic
models. `scripts/validate_schemas.py` fails when a checked-in schema differs from the current
model.

| Schema | Version | Unit |
|---|---|---|
| Candidate evidence bundle | 1.1 | One query-candidate record per JSONL line |
| Run manifest | 1.1 | One complete run provenance record |
| Structure prediction queue | 1.0 | One manually reviewed pair proposal |

Version 1.1 adds optional candidate identity/policy metadata and normalization provenance, so
1.0 readers that tolerate unknown optional fields remain compatible. The candidate TSV is a
derived flat view with a fixed column order in `outputs/candidates.py`.

## Score semantics

Scores are optional floats from 0.0 to 1.0. `None` means not calculated or unavailable.
`0.0` means a calculation was performed and produced zero. Total score is not calculated by
model validation, and an Evidence Tier is not inferred from a score.

## Evidence semantics

`EvidenceStatus.not_run` represents every MVP-1A engine that has not run. `missing` means the
engine or input was applicable but no usable evidence was available. Exact protein evidence
and `ortholog_transferred` evidence are represented by different `EvidenceOrigin` values and
must not be merged silently. Provenance and warnings remain attached to each evidence record.

## Compatibility

- Patch: documentation corrections or optional descriptive metadata that does not affect
  existing consumers.
- Minor: backward-compatible optional fields or enum additions documented for tolerant readers.
- Major: removed/renamed required fields, changed meaning/type, or other breaking changes.

No migration engine is included in MVP-1A.

## Identifier and candidate semantics

Original and normalized identifiers coexist in each candidate snapshot. `exact_match`,
`unique_alias_match`, `ambiguous_match`, and `no_match` are separate states. Missing coordinate
and missing annotation are independent booleans/warnings and are not encoded as empty numeric
values. Candidate list fields use sorted `|`-joined values in TSV and arrays in canonical
JSONL.
