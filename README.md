# ProteinInteractionHunter

ProteinInteractionHunter is a local-first research support system intended to prioritize
proteins that may physically interact with, or be functionally associated with, one or more
query proteins in the same organism. It does **not** establish that a physical interaction
exists and does not replace experimental validation.

This repository is independent from ProteinHunter. It does not import, depend on, copy source
from, or modify ProteinHunter_v5.

## Current status: MVP-1

MVP-1 resolves query identifiers, enumerates every proteome record for each query, and evaluates
local gene context, operon proxy, domain, functional-complementarity, localization, orthology,
phylogenetic-profile, gene-fusion, and known-interaction evidence. Optional integrated scoring,
per-query dense ranking, and evidence tier classification remain fully auditable through raw
evidence and component/category breakdowns.

Scores, ranks, and Evidence Tiers prioritize candidates; they do **not** establish a physical
interaction and do not replace experimental validation. The machine-readable source of truth is
the UTF-8 JSONL candidate evidence bundle. TSV and Excel are derived views.

## Scientific and data policy

- Physical interaction is never asserted from coordinate proximity.
- AlphaFold 3, AlphaFold Server, ColabFold, and other structure predictors are not run or
  submitted to automatically.
- No external service or network client is enabled in MVP-1; execution is local-only.
- Private sequences are not transmitted externally.
- Hypothetical proteins are not excluded solely because they lack annotation.
- Absence of evidence is not treated as negative evidence.
- `None` means a score has not been calculated; `0.0` is a calculated zero.

## Installation

Python 3.12 or newer is required.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

Development installation:

```bash
python -m pip install -e '.[dev]'
```

## CLI

```bash
protein-interaction-hunter --version
protein-interaction-hunter validate-config --config tests/fixtures/config.valid.yaml
protein-interaction-hunter validate-inputs --config tests/fixtures/config.valid.yaml
protein-interaction-hunter inspect-fixture --config tests/fixtures/config.valid.yaml
protein-interaction-hunter generate-candidates --config tests/fixtures/config.valid.yaml
```

Input paths in YAML are resolved relative to the YAML file. Evidence engines are optional and
default to disabled when their optional section is omitted. `gene_context.enabled` controls both
gene-context and operon-proxy evaluation; `scoring.enabled` controls integrated scoring, and
`evidence_tiers.enabled` requires scoring. `require_query_coordinates: true` makes missing or
ambiguous query coordinates fatal.

`generate-candidates` writes `candidate_evidence_bundle.jsonl`, `candidate_table.tsv`,
`ProteinInteractionHunter.xlsx`, `run_manifest.json`, `config.snapshot.yaml`, and
`warning_summary.tsv` under `output.directory`. When enabled, scoring and tier results are
included in JSONL and TSV and in the `Integrated_Scoring`, `Scoring_Components`,
`Candidate_Ranking`, `Evidence_Tiers`, and `Tier_Summary` Excel sheets.

## Coordinate and context rules

GFF coordinates remain 1-based closed intervals. The parser retains seqid, source, feature
type, start/end, strand, IDs, all parents, attributes, warnings, and optional
`##sequence-region` boundaries. Unknown strand is `?`; reversed coordinates are invalid.
Circular origin wrapping and GFF phase interpretation are not implemented.

For each protein, a unique CDS coordinate is the deterministic representative. A related wider
gene extent is retained as a warning; multiple distinct CDS coordinates are ambiguous rather
than selected by row order. Gene/CDS rows for one biological unit collapse before feature
counting. Unmapped gene/RNA units remain visible to all-feature counts.

For same-contig pairs, separated distance is `max(start) - min(end) - 1`; overlap has distance
zero and inclusive width `min(end) - max(start) + 1`. Different-contig distance/overlap are null
with `not_applicable`. Missing and ambiguous coordinates remain null with `missing` and `failed`.
Genomic left/right is separate from upstream/downstream in the query's transcription direction.
Opposite non-overlapping strands are convergent (`+ … -`) or divergent (`- … +`); overlapping
opposite strands are `opposite_parallel`.

Neighborhood membership is based only on absolute representative feature-index delta. Contig
edge distances require `##sequence-region`; otherwise they and completeness are unknown.
The rule/provenance version is `mvp1b-gene-context-v1`.

## Identifier normalization

FASTA protein IDs remain canonical and their original spellings are never replaced. Matching
aliases are NFKC-normalized, trimmed, case-folded, and stripped of a leading known namespace
prefix (`protein[_id]`, `gene[_id]`, `locus[_tag]`, `old_locus_tag`, `cds`, `id`, or `parent`)
when followed by `:` or `=`. Protein ID, gene ID, locus tag, old locus tag, GFF `ID`/`Parent`,
and `protein_id` participate in an auditable many-to-many index. Exact, unique-alias,
ambiguous, and absent matches are distinct; ambiguous query IDs are fatal.

## FASTA fixture header convention

The first whitespace-delimited token is the protein ID. Free text is the description.
Optional metadata uses bracketed keys:

```text
>QUERY_001 Query enzyme [gene_id=gene_query] [locus_tag=LT0001]
```

Supported optional keys are `gene_id` and `locus_tag`. Duplicate IDs are invalid. Duplicate
sequences remain separate loci and receive a duplicate group.

## Repository layout

```text
src/protein_interaction_hunter/  package, domain models, services, adapters, writers
schemas/                         checked-in Draft 2020-12 JSON Schemas
tests/fixtures/                  artificial local-only validation data
docs/                            design, architecture, schemas, development
scripts/                         schema drift validation
```

See the [detailed design](docs/ProteinInteractionHunter_DESIGN.md),
[architecture](docs/architecture.md), [schemas](docs/schemas.md), and
[development guide](docs/development.md).

## Limitations and roadmap

MVP-1 remains a local evidence-prioritization workflow. It does not submit structures,
retrieve network evidence, infer circular-origin wrapping, or replace manual review and
experimental validation. Future work may add isolated external adapters and structure-result
import without changing the meaning of current raw evidence.
