# Development

## Environment

Use Python 3.12 on WSL/Linux:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
```

## Checks

```bash
python -m pytest
python -m ruff check .
python -m mypy src
python scripts/validate_schemas.py
python -m pytest --cov=protein_interaction_hunter --cov-report=term-missing
```

Focused MVP-1 integration checks:

```bash
python -m pytest tests/integration/test_candidate_pipeline.py
python -m pytest tests/integration/test_final_integration_audit.py
protein-interaction-hunter generate-candidates --config tests/fixtures/config.all_engines.yaml
```

The gene-context tests must cover adjacency, inclusive overlap, self pairs, transcription-relative
position, all strand relations, representative feature de-duplication, ambiguous coordinates,
different contigs, missing coordinates, sequence-region boundaries, deterministic output, and
disabled/required-coordinate policy behavior. Update checked-in JSON schemas with
`python scripts/validate_schemas.py --write` only after an intentional model/version change.

Before committing, also run:

```bash
git status --short
git diff --check
git diff --stat
```

The repository is stored on a Windows-mounted path. Keep text files UTF-8 with LF endings,
avoid unnecessary executable bits, and use repository-local `core.autocrlf=false` and
`core.filemode=false`; do not change global Git configuration.
