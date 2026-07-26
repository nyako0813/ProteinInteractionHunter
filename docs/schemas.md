# Schemas

All checked-in JSON Schemas use JSON Schema Draft 2020-12 and are generated from the Pydantic
models. `scripts/validate_schemas.py` fails when a checked-in schema differs from the current
model.

| Schema | Version | Unit |
|---|---|---|
| Candidate evidence bundle | 1.3 | One query-candidate record per JSONL line |
| Run manifest | 1.3 | One complete run provenance record |
| Structure prediction queue | 1.0 | One manually reviewed pair proposal |

Version 1.3 adds the MVP-1L evidence-tier result and manifest configuration provenance.
These are additive changes; the legacy EvidenceTier enum values remain parseable for backward
compatibility. The candidate TSV is a derived flat view with a fixed column order in
`outputs/candidates.py`.

## Gene-context semantics

Coordinates are 1-based closed intervals. For same-contig intervals:

- separated: `distance_bp = max(start) - min(end) - 1`, `overlap_bp = 0`;
- overlapping: `distance_bp = 0`, `overlap_bp = min(end) - max(start) + 1`;
- same feature: explicit `same_feature`, with the inclusive self-overlap retained;
- different contigs: both numeric values are null and status is `not_applicable`;
- missing/ambiguous coordinates: numeric values are null and status is `missing`/`failed`.

`coordinate_position` is genomic left/right/overlap and does not depend on strand.
`relative_position` is upstream/downstream in the query's transcription direction. Unknown
query strand yields `unknown`. `strand_relationship` is `same_direction`, `convergent`
(left `+`, right `-`), `divergent` (left `-`, right `+`), `opposite_parallel` for overlapping
opposite-strand intervals, `unknown`, or `different_contig`.

Intervening counts exclude the query and candidate. For separated pairs, a representative unit
is counted when it lies wholly or partly inside the open interval between the pair; overlap and
same-feature pairs have zero intervening units. `intervening_gene_count` counts gene/CDS units,
while `intervening_feature_count` also retains other representative units such as RNA.

`feature_index_delta` is the absolute difference between deterministic representative indices.
`within_neighborhood_window` is true when that delta is no greater than
`gene_context.neighborhood_gene_count`. Completeness reports whether the configured index window
can be observed on both sides: `complete`, `left_truncated`, `right_truncated`,
`both_truncated`, or `unknown` when contig boundaries are unavailable.

## Evidence and score semantics

`available` means coordinate facts were calculated. `missing` means an applicable coordinate
was unavailable. `failed` includes ambiguous coordinate mapping. `not_applicable` is used for
different contigs, and `not_run` for disabled/unimplemented engines. Missing evidence is not
negative evidence.

Backward-compatible candidate summary scores are optional floats from 0.0 to 1.0. `None`
means not calculated or unavailable; `0.0` is a calculated zero. The integrated scoring record
also stores its configured output scale, component/category breakdown, sufficiency, and rank.
Evidence tiers are optional prioritization labels and are not interaction claims.

## Compatibility

- Patch: documentation corrections or optional descriptive metadata with no consumer impact.
- Minor: backward-compatible optional fields or enum additions documented for tolerant readers.
- Major: removed/renamed required fields, changed meaning/type, or other breaking changes.

No automatic schema migration engine is included.
