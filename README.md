# ProteinInteractionHunter

ProteinInteractionHunter is a local-first research support system intended to prioritize
proteins that may physically interact with, or be functionally associated with, one or more
query proteins in the same organism. It does **not** establish that a physical interaction
exists and does not replace experimental validation.

This repository is independent from ProteinHunter. It does not import, depend on, or copy
source code from ProteinHunter_v5.

## Current status: MVP-0

MVP-0 provides project structure, strict configuration validation, domain models, versioned
schemas, synthetic fixtures, local parsers, manifest generation, and schema-only output
writers. Candidate generation, ranking, scoring, evidence tiers, BLAST/DIAMOND, annotation
retrieval, and biological analysis engines are intentionally not implemented. Calling the
pipeline raises `NotImplementedError`; it never returns a fabricated successful analysis.

The machine-readable source of truth is the UTF-8 JSONL candidate evidence bundle. TSV and
Excel reports are derived views and never replace the JSONL record.

## Scientific and data policy

- Physical interaction is never asserted solely from a score or tier.
- AlphaFold 3, AlphaFold Server, ColabFold, and other structure predictors are not run or
  submitted to automatically.
- No external service or network client is included in MVP-0; validation is local-only.
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
```

The input paths in a YAML file are resolved relative to the directory containing that YAML
file, not relative to the current shell directory. Config-only validation does not require
optional files to exist; `validate-inputs` checks required files and formats.

## FASTA fixture header convention

The first whitespace-delimited token is the protein ID. Free text is the description.
Optional metadata uses bracketed keys:

```text
>QUERY_001 Query enzyme [gene_id=gene_query] [locus_tag=LT0001]
```

Supported optional keys are `gene_id` and `locus_tag`. Duplicate IDs are invalid. Duplicate
sequences are valid and reported as groups so a later MVP can flag them without discarding
their distinct loci.

## Repository layout

```text
src/protein_interaction_hunter/  package, domain models, ports, adapters, writers
schemas/                         checked-in Draft 2020-12 JSON Schemas
tests/fixtures/                  artificial local-only validation data
docs/                            design, architecture, schemas, development
scripts/                         schema drift validation
```

See the [detailed design](docs/ProteinInteractionHunter_DESIGN.md),
[architecture](docs/architecture.md), [schemas](docs/schemas.md), and
[development guide](docs/development.md).

## Limitations and roadmap

MVP-1 will add local single-organism candidate generation, identifier normalization, genomic
distance and context evidence, transparent scoring traces, contradictions, derived reports,
and a manual structure-prediction queue. MVP-2 will add optional evolutionary evidence.
MVP-3 will isolate optional external evidence behind adapters. Any future structure-result
import will remain separate from the initial non-structural ranking.
