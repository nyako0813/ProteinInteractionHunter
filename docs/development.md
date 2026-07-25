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

Before committing, also run:

```bash
git status --short
git diff --check
git diff --stat
```

The repository is stored on a Windows-mounted path. Keep text files UTF-8 with LF endings,
avoid unnecessary executable bits, and use repository-local `core.autocrlf=false` and
`core.filemode=false`; do not change global Git configuration.
