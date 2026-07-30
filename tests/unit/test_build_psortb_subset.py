import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from scripts.build_psortb_subset import build_subset, write_subset  # noqa: E402


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="\n")
    return path


def test_subset_is_deterministic_reasoned_and_accession_only(tmp_path: Path) -> None:
    fasta = _write(
        tmp_path / "proteome.faa",
        ">Q query description\nMMMMMMMMMMMMMMMMMMMM\n"
        ">A membrane product\nLLLLLLLLLLLLLLLLLLLL\n"
        ">B no hit\nDEKRDEKRDEKRDEKRDEKR\n",
    )
    seed = _write(
        tmp_path / "seed.tsv",
        "protein_id\tlength_aa\tdescription\tselection_reasons\n"
        "Q\t20\tquery description\tquery|query_flank_10_cds\n",
    )
    domains = _write(
        tmp_path / "domains.tsv",
        "protein_id\tsource\taccession\tname\tstart\tend\tarchitecture_index\n"
        "A\tPfam\tPF1\tone\t1\t5\t0\n"
        "A\tPfam\tPF2\ttwo\t10\t15\t1\n",
    )

    proteins, reasons = build_subset(
        proteome_fasta=fasta,
        seed_audit=seed,
        domain_table=domains,
        sample_size=1,
    )
    output = tmp_path / "subset.faa"
    audit = tmp_path / "audit.tsv"
    write_subset(proteins, reasons, fasta_output=output, audit_output=audit)

    assert [protein.protein_id for protein in proteins] == ["A", "B", "Q"]
    assert "multiple_pfam_domains" in reasons["A"]
    assert "no_pfam_hit" in reasons["B"]
    assert output.read_text(encoding="utf-8").splitlines()[::2] == [">A", ">B", ">Q"]
    assert "description" not in output.read_text(encoding="utf-8")
    with audit.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["protein_id"] for row in rows] == ["A", "B", "Q"]
    assert b"\r" not in output.read_bytes()
    assert b"\r" not in audit.read_bytes()
