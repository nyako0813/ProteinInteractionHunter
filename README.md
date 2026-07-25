# ProteinInteractionHunter

ProteinInteractionHunter is a local-first research support system intended to prioritize
proteins that may physically interact with, or be functionally associated with, one or more
query proteins in the same organism. It does **not** establish that a physical interaction
exists and does not replace experimental validation.

This repository is independent from ProteinHunter. It does not import, depend on, copy source
from, or modify ProteinHunter_v5.

## Current status: MVP-1B

MVP-1B resolves one or more query identifiers, enumerates every proteome record for each query,
normalizes GFF3 coordinates, and writes observable gene-context facts for every query-candidate
pair. It reports same/different contig, interval distance and overlap, coordinate and
transcription-relative position, strand relationship, intervening representative features,
feature-index neighborhood membership, contig-edge distances, completeness, status, warnings,
and provenance.

This release does **not** infer operons, conserved neighborhoods, orthologs, functional
relationships, physical interactions, ranks, scores, Evidence Tiers, or structure results.
The machine-readable source of truth is the UTF-8 JSONL candidate evidence bundle. TSV and
Excel are derived views and never replace the JSONL record.

## Scientific and data policy

- Physical interaction is never asserted from coordinate proximity.
- AlphaFold 3, AlphaFold Server, ColabFold, and other structure predictors are not run or
  submitted to automatically.
- No external service or network client is included in MVP-1B; execution is local-only.
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

Input paths in YAML are resolved relative to the YAML file. `gene_context.enabled` controls the
MVP-1B engine, `neighborhood_gene_count` controls the representative-index window, and
`require_query_coordinates: true` makes missing or ambiguous query coordinates fatal.

`generate-candidates` writes `candidate_evidence_bundle.jsonl`, `candidate_table.tsv`,
`ProteinInteractionHunter.xlsx`, `run_manifest.json`, `config.snapshot.yaml`, and
`warning_summary.tsv` under `output.directory`. It does not calculate or print ranks, scores,
top candidates, relationship claims, or Evidence Tiers.

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

Later MVP-1 stages will add operon proxies, functional rules, transparent scoring traces,
contradictions, and a manual structure-prediction queue. MVP-2 will add optional orthology,
phylogenetic-profile, conserved-neighborhood, and fusion evidence. MVP-3 will isolate optional
external evidence behind adapters. Structure-result import will remain separate from the
initial non-structural ranking.
