# MVP-1 Final Integration Audit

## Decision

**READY_WITH_LIMITATIONS**

This decision means the local MVP-1 pipeline is ready for a small, supervised real-data pilot.
It is not a production-readiness statement and does not change the scientific meaning of scores,
ranks, or Evidence Tiers.

## Baseline and environment

- Baseline commit: `3da7e45d9e2698f3f8dd974978baf0b3ba20991a`
  (`Implement MVP-1L evidence tier classification`)
- Audit date: 2026-07-26
- Environment: WSL Ubuntu, Python 3.12, branch `main`
- Start state: clean; `HEAD == origin/main`
- Audit fixture: `tests/fixtures/config.all_engines.yaml`
- Canonical fixture size: 1 query, 13 proteins, 13 query-candidate pairs
- `_UNIMPLEMENTED_ENGINES`: empty tuple

The audit covers candidate generation, gene context, operon proxy, domain rules, functional
complementarity, annotation-only localization, orthology, phylogenetic profiles, gene fusion,
known interactions, integrated scoring, per-query ranking, evidence tiers, JSONL, TSV, Excel,
manifest, warnings, schemas, config validation, and CLI behavior.

## Execution matrix

| Mode | Outcome |
|---|---|
| All engines enabled | 13 bundles; all enabled engine statuses were non-`not_run` |
| All optional engines disabled | Candidate generation completed; raw evidence/scoring/tier empty; statuses `not_run` |
| Raw evidence enabled, scoring disabled | Raw evidence identical to canonical; no integrated score/rank/tier |
| Raw evidence + scoring enabled, tiers disabled | Raw evidence, score, rank, components, and categories identical; tier null/result empty |
| Section omitted | Sections with backward-compatible defaults (profile, fusion, known interactions, scoring, tiers) behave as disabled; required core sections fail config validation when omitted |
| Enabled table/rules path absent | Fail-fast `InputValidationError` or loader validation error |
| Malformed required row/column/value | Fail-fast with a field-specific validation message |
| Header-only domain/orthology/profile table | Per-pair `missing`; scoring component excluded from denominator |
| Header-only fusion/known-interaction table | Evaluated `available` record-absence; neutral component value 0 |
| Partially or ambiguously mapped record | Missing/uncertain or positive attenuation; never promoted to strong positive |

The fail-fast policy is consistent for malformed local input. The pipeline does not currently
isolate a malformed optional engine and continue other engines.

## Canonical all-engine E2E

The pipeline completed with 13 deterministic candidate rows and one bundle per pair. The self
pair remains present, is `excluded`, has no rank, and is `unclassified` with `self_pair` in
`failed_requirements`. Every pair has a scoring result and tier result when both engines are
enabled. Eligible, sufficient, non-excluded pairs have a formal score and rank.

Default tier distribution:

| Tier | Count |
|---|---:|
| Tier 1 | 1 |
| Tier 2 | 0 |
| Tier 3 | 2 |
| Tier 4 | 9 |
| Unclassified | 1 |

## Determinism and input order

Three runs with the same config and inputs were compared. Candidate TSV, evidence JSONL, warning
summary, every Excel sheet value and row order, candidate order, scores, ranks, tiers,
component/category order, terms, warnings, and provenance were identical. Only manifest
timestamps were normalized before comparison.

Reversing rows independently in annotation, domain, orthology, phylogenetic-profile, fusion, and
known-interaction tables produced identical pair semantics, raw evidence, score, rank, tier,
terms, warnings, and candidate order after normalizing the config-derived run ID. FASTA was not
reordered because FASTA order is the canonical protein enumeration input.

## Individual engine ablation

Each configurable evidence engine was disabled separately:

- gene context (and its coupled operon proxy);
- domain;
- functional complementarity;
- localization;
- orthology;
- phylogenetic profile;
- fusion;
- known interactions.

Only the selected raw evidence fields changed. Candidate order, candidate content, disposition,
and relationship remained identical. Disabled score components were present with `applied=false`,
`evidence_status=not_run`, and `exclusion_reason=evidence_not_run`. Effective-weight and
contribution recomputation confirmed they were removed from the denominator rather than treated
as zero. Resulting score/rank/tier changes were explainable by that component/category or
high-specificity removal.

Operon does not have an independent enabled flag: it is intentionally coupled to
`gene_context.enabled`. Therefore an operon-only ablation is not available in MVP-1.

## Duplicate resistance and pair symmetry

Exact duplicates in domain, orthology, phylogenetic-profile, fusion, and known-interaction local
tables are rejected explicitly by their loaders. Swapped duplicates are also rejected where the
pair is canonicalized. Duplicate records therefore cannot double a score, high-specificity
count, publication count, source count, or support term.

Known-interaction and fusion canonical-pair symmetry tests passed. Orthology support is symmetric
while retaining query/candidate-oriented output fields. The audit found that domain rules were
evaluated only in the declared query-role direction. The evaluator now tests both role
orientations while retaining oriented protein/accession fields; forward and reverse pairs have
the same match/rule semantics.

## Identifier mapping

Existing and audit tests cover exact protein ID, gene/locus aliases, normalization of case,
whitespace and known prefixes, unknown identifiers, ambiguous aliases, duplicate canonical
query resolution, and reversed local pair orientation.

- Exact and unique aliases preserve canonical FASTA protein IDs.
- Ambiguous query resolution is fatal rather than guessed.
- Candidate identifier ambiguity attenuates applied positive scoring components.
- Uncertain known-interaction mapping does not satisfy direct/physical high-specificity gates.
- Unknown or missing mapping is unavailable evidence, not negative evidence.
- Mapping warnings remain in bundle and warning-summary output.

## Missing values and numeric boundaries

Loader and model tests cover blank optional values, missing annotation fields, missing
coordinates, missing domain/orthology/profile records, blank fusion numeric fields, missing known
interaction method/confidence, and missing source metadata. Canonical TSV/JSON/YAML outputs
contained no inappropriate literal `nan`, `None`, or `null` cells and no CRLF. Excel empty values
were cells; the literal `none` observed in localization output is the fixture's valid
transmembrane annotation value.

Boundary tests cover inclusive coordinate overlap/distance, exact neighborhood threshold,
different contigs, reversed coordinates, phylogenetic similarity and informative/discordance
conditions, fusion coverage/overlap/support thresholds, malformed nonfinite scoring input,
category caps, zero available weight, exact sufficiency boundaries, exact tier thresholds,
negative limits, high-specificity counts, and tier caps.

## Independent recomputation

All 13 canonical pairs passed independent score recomputation:

- each contribution equals `normalized_value * effective_weight`;
- applied contributions sum to `raw_weighted_sum`;
- effective weights sum to `available_weight`;
- clamped normalized and provisional values match at 12-decimal calculation precision;
- formal score presence matches sufficiency.

All ranks were independently recomputed per query using formal score descending, configured tie
precision, dense ranks, and deterministic secondary display keys. Excluded and insufficient
pairs have no rank.

All tiers were independently recomputed from formal score, sufficiency, category/component
counts, weight, negative count, raw direct/physical/fusion flags, predicted-only,
functional-association-only, base-tier gates, and caps. Base tier, assigned tier, caps, and
Unclassified reason matched for every pair.

## Cross-output consistency

Internal models, schema-validated JSONL, 213-column Candidate TSV, Candidate Ranking,
Integrated Scoring, Scoring Components, Evidence Tiers, and Tier Summary were compared by pair.
IDs, disposition, relationship, scores, ranks, tiers, sufficiency, counts, terms, warnings, and
rule versions matched. The workbook contained 13 candidate/scoring/tier rows, 117 scoring
component rows, and one query tier-summary row.

## Schema and manifest

All three schemas were regenerated from Pydantic models and checked for zero drift. Draft
2020-12 validation succeeded for checked-in schemas and canonical JSONL/manifest output. Bundle
and manifest schema versions are 1.3. Legacy A-X EvidenceTier enum values remain parseable.
Unknown model fields remain rejected.

Input SHA256 and size were stable on repeated reads; changing one character changed only the
target input hash. Generated outputs are not registered as inputs. The manifest records
software/schema/rule versions, input hashes and sizes, config snapshot path, scoring config, tier
config, warnings, and incomplete-evidence flags.

Limitation: per-pair engine statuses live in the canonical bundle, and raw-engine enabled state
lives in the full config snapshot; the run manifest does not duplicate a run-level engine-status
matrix or an output-artifact checksum list.

## Warning audit

Canonical warning-summary counts exactly matched the multiset of bundle warnings, while manifest
warnings matched the unique warning codes. Disabled engines do not generate errors. Routine
`unclassified` tier outcomes now retain their reason in `failed_requirements` without creating a
warning. Missing optional biological coverage remains explicit but is not converted to negative
evidence.

## CLI audit

The following real entry-point commands succeeded:

```text
python -m protein_interaction_hunter --help
python -m protein_interaction_hunter --version
python -m protein_interaction_hunter validate-config --config <config>
python -m protein_interaction_hunter validate-inputs --config <config>
python -m protein_interaction_hunter generate-candidates --config <config>
```

The performance fixture reported 150 proteins, 150 pairs, 149 included, and one excluded self
pair consistently with generated output. Invalid config tests return nonzero with a concise
field-specific error and no traceback. Obsolete MVP-0/MVP-1B and “tiers remain unset” help text
was removed.

## Performance smoke

One query against 150 synthetic proteins was run with all local engines enabled:

| Metric | Result |
|---|---:|
| Candidate pairs | 150 |
| Runtime | 1.534 s |
| Peak RSS | 73,628 KiB |
| Total output | 3,884,302 bytes |
| Workbook candidate/scoring/tier rows | 150 each |
| Scoring component rows | 1,350 |

No extreme runtime, memory growth, loop, or workbook expansion was observed. This is a smoke
test, not a benchmark or scaling guarantee.

## Final validation gates

- `pytest -q`: **324 passed** in 63.79 seconds.
- Final-integration plus domain/functional symmetry tests: **45 passed**.
- `ruff check src tests`: passed.
- `mypy src tests`: passed with no issues in 107 source files.
- `scripts/validate_schemas.py --write`: regenerated three schemas.
- `scripts/validate_schemas.py`: validated three schemas with no drift.
- `git diff --check`: passed.
- Real CLI and 150-protein performance smoke: passed.

## Scientific interpretation

### NEAR_001

- Formal score `80.6896551724`, rank 1, base/assigned Tier 1.
- All nine components are positive.
- Direct known-interaction and fusion-association evidence provide two high-specificity classes.
- Direct interaction category weight is capped; no negative component or tier cap is present.
- This is strong software-internal evidence, not a confirmed interaction.

### MID_001

- Formal score `16.6666666667`, rank 4, base/assigned Tier 4.
- Known functional association is the only positive component.
- Direct, physical, fusion, and high-specificity support are absent.
- `functional_association_only=true`; the default Tier 3 cap is non-operative because Tier 4 is
  already lower. A controlled cap-boundary integration case lowers the quantitative gates and
  confirms base Tier 2 is capped to assigned Tier 3.

### FRAG_001

- Formal score `5.7142857143`, rank 11, base/assigned Tier 4.
- The known-interaction record is predicted-only and supplies no high-specificity class.
- `predicted_only=true`; the default Tier 3 cap is non-operative because Tier 4 is already lower.
  A controlled cap-boundary test confirms base Tier 2 is capped to Tier 3.

### Self pair

- The candidate row is retained with `excluded` disposition.
- Formal score is auditable, but rank is null.
- Base and assigned tier are `unclassified`; `failed_requirements=["self_pair"]`.
- The routine Unclassified result is not a warning.

## Defects found and fixed

1. **Routine Unclassified outcomes were always warnings.** Cause: eligibility and quantitative
   gate results populated both `failed_requirements` and `warnings`. Fix: keep classification
   reasons in `failed_requirements`; reserve warnings for evaluation anomalies.
2. **Domain-pair evaluation was direction-dependent.** Cause: rules were tested only as
   query-role to candidate-role. Fix: evaluate both role orientations while retaining oriented
   output fields.
3. **User-facing capability documentation was stale.** Cause: CLI, README, architecture,
   schema, development text, and a no-config error still described MVP-0/MVP-1B. Fix: update only
   affected capability and boundary text; scientific caution remains explicit.
4. **The required whole-suite Mypy gate failed.** Cause: several test helpers lacked return or
   parameter annotations, and the new audit used untyped deserialization results. Fix: add
   typing-only annotations/casts; no runtime behavior changed.

No score weights, category caps, tier thresholds, fixture expectations, or raw evidence meanings
were changed.

## Remaining limitations

- Localization is annotation-only; there is no independent localization-table adapter.
- Operon proxy enablement is coupled to gene context.
- Malformed optional engine input is fail-fast rather than isolated.
- Manifest does not embed per-engine aggregate statuses or output artifact hashes.
- The synthetic fixture has one query; rank partitioning across multiple queries is covered by
  scoring unit tests rather than this canonical E2E.
- No network retrieval, conserved-neighborhood engine, circular-origin inference, or automatic
  structure submission is included.

Within these limits, the project is ready for a small, supervised real-data pilot.
