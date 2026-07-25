"""Reasonable static checks for repository separation and no network clients."""

import ast
from pathlib import Path

FORBIDDEN_IMPORT_ROOTS = {
    "Bio",
    "ProteinHunter",
    "alphafold",
    "biogrid",
    "colabfold",
    "httpx",
    "intact",
    "requests",
    "stringdb",
}


def test_package_has_no_forbidden_runtime_imports_or_external_path() -> None:
    root = Path(__file__).resolve().parents[2]
    source_root = root / "src" / "protein_interaction_hunter"
    imported_roots: set[str] = set()
    for path in source_root.rglob("*.py"):
        text = path.read_text(encoding="utf-8")
        assert "/home/nyako/projects/ProteinHunter_v5" not in text
        tree = ast.parse(text, filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported_roots.update(alias.name.split(".", 1)[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported_roots.add(node.module.split(".", 1)[0])
    assert imported_roots.isdisjoint(FORBIDDEN_IMPORT_ROOTS)


def test_project_dependencies_are_local_only() -> None:
    root = Path(__file__).resolve().parents[2]
    pyproject = (root / "pyproject.toml").read_text(encoding="utf-8").lower()
    for dependency in ("requests", "httpx", "biopython", "blast", "diamond", "colabfold"):
        assert f'"{dependency}' not in pyproject
