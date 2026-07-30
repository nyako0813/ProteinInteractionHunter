"""Run an auditable per-species DIAMOND homolog and reciprocal-hit pilot."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

OUTPUT_COLUMNS = (
    "query_species",
    "target_species",
    "query_protein_id",
    "target_protein_id",
    "percent_identity",
    "alignment_length",
    "query_length",
    "target_length",
    "query_coverage",
    "target_coverage",
    "evalue",
    "bit_score",
    "rank",
    "best_hit_tie",
    "reciprocal_best_hit",
    "reciprocal_tie",
    "length_ratio",
    "fragment_flag",
    "paralog_count",
    "decision",
    "decision_reason",
    "rule_version",
)


@dataclass(frozen=True)
class Hit:
    query_id: str
    subject_id: str
    percent_identity: float
    alignment_length: int
    query_length: int
    subject_length: int
    evalue: float
    bit_score: float

    @property
    def query_coverage(self) -> float:
        return self.alignment_length / self.query_length

    @property
    def subject_coverage(self) -> float:
        return self.alignment_length / self.subject_length

    @property
    def length_ratio(self) -> float:
        return self.subject_length / self.query_length


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    current: str | None = None
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if current is not None:
                records[current] = "".join(chunks)
            current = line[1:].split()[0]
            chunks = []
        elif line.strip():
            chunks.append(line.strip())
    if current is not None:
        records[current] = "".join(chunks)
    return records


def _parse_hits(path: Path) -> list[Hit]:
    hits: list[Hit] = []
    if not path.exists():
        return hits
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) != 8:
            raise ValueError(f"Malformed DIAMOND row {path}:{line_number}")
        hit = Hit(
            query_id=parts[0],
            subject_id=parts[1],
            percent_identity=float(parts[2]),
            alignment_length=int(parts[3]),
            query_length=int(parts[4]),
            subject_length=int(parts[5]),
            evalue=float(parts[6]),
            bit_score=float(parts[7]),
        )
        if not all(
            math.isfinite(value) for value in (hit.percent_identity, hit.evalue, hit.bit_score)
        ):
            raise ValueError(f"Non-finite DIAMOND value {path}:{line_number}")
        hits.append(hit)
    return sorted(hits, key=lambda hit: (-hit.bit_score, hit.evalue, hit.subject_id))


def _passes(hit: Hit, search: dict[str, Any]) -> bool:
    return (
        hit.evalue <= float(search["evalue_maximum"])
        and hit.bit_score >= float(search["minimum_bit_score"])
        and hit.percent_identity >= float(search["minimum_percent_identity"])
        and hit.alignment_length >= int(search["minimum_alignment_length"])
        and hit.query_coverage >= float(search["minimum_query_coverage"])
        and hit.subject_coverage >= float(search["minimum_subject_coverage"])
        and float(search["minimum_length_ratio"])
        <= hit.length_ratio
        <= float(search["maximum_length_ratio"])
    )


def classify_species(
    hits: list[Hit],
    reciprocal_hits: dict[str, list[Hit]],
    *,
    query_normalized_id: str,
    policy: dict[str, Any],
) -> tuple[str, str, bool, bool, int]:
    search = policy["search"]
    classification = policy["classification"]
    if not hits:
        return "no_detectable_homolog", "no_hit_at_search_evalue", False, False, 0
    passing = [hit for hit in hits if _passes(hit, search)]
    top = hits[0]
    fragment = (
        top.length_ratio <= float(classification["fragment_length_ratio_maximum"])
        or top.query_coverage < float(search["minimum_query_coverage"])
        or top.subject_coverage < float(search["minimum_subject_coverage"])
    )
    if not passing:
        if fragment:
            return "fragment_only", "only_fragmentary_or_short_coverage_hits", False, False, 0
        return "weak_homolog", "hits_fail_conservative_candidate_thresholds", False, False, 0
    best = passing[0]
    tie_fraction = float(search["best_hit_tie_relative_bit_score"])
    best_tie = sum(hit.bit_score >= best.bit_score * (1.0 - tie_fraction) for hit in passing) > 1
    secondary_fraction = float(search["secondary_hit_relative_bit_score"])
    paralog_count = max(
        0,
        sum(hit.bit_score >= best.bit_score * secondary_fraction for hit in passing) - 1,
    )
    reciprocal = reciprocal_hits.get(best.subject_id, [])
    reciprocal_passing = [hit for hit in reciprocal if _passes(hit, search)]
    reciprocal_best = bool(
        reciprocal_passing and reciprocal_passing[0].subject_id == query_normalized_id
    )
    reciprocal_tie = False
    if reciprocal_passing:
        reciprocal_top = reciprocal_passing[0].bit_score
        reciprocal_tie = (
            sum(
                hit.bit_score >= reciprocal_top * (1.0 - tie_fraction) for hit in reciprocal_passing
            )
            > 1
        )
    if best_tie or reciprocal_tie:
        return (
            "ambiguous_RBH",
            "forward_or_reciprocal_best_hit_tie",
            reciprocal_best,
            reciprocal_tie,
            paralog_count,
        )
    if paralog_count >= int(classification["multi_copy_minimum_secondary_hits"]):
        return (
            "multi_copy_homolog_family",
            "near_top_secondary_candidates_present",
            reciprocal_best,
            reciprocal_tie,
            paralog_count,
        )
    if reciprocal_best:
        return (
            "unique_RBH",
            "conservative_thresholds_and_unique_reciprocal_best_hit",
            True,
            False,
            paralog_count,
        )
    return (
        "candidate_homolog",
        "conservative_thresholds_met_without_unique_reciprocal_support",
        False,
        reciprocal_tie,
        paralog_count,
    )


def _run(command: list[str], stdout_path: Path | None = None) -> None:
    if stdout_path is None:
        result = subprocess.run(command, check=False, capture_output=True, text=True)
    else:
        with stdout_path.open("w", encoding="utf-8", newline="\n") as handle:
            result = subprocess.run(
                command, check=False, stdout=handle, stderr=subprocess.PIPE, text=True
            )
    if result.returncode:
        raise RuntimeError(
            f"Command failed ({result.returncode}): {' '.join(command)}\n{result.stderr}"
        )


def run_pilot(
    *,
    diamond_path: Path,
    normalized_directory: Path,
    panel_path: Path,
    mapping_path: Path,
    policy_path: Path,
    working_directory: Path,
    output_path: Path,
    summary_path: Path,
    metadata_path: Path,
) -> dict[str, Any]:
    policy = yaml.safe_load(policy_path.read_text(encoding="utf-8"))
    query_raw_id = policy["query_protein_id"]
    with panel_path.open(encoding="utf-8", newline="") as handle:
        panel = list(csv.DictReader(handle, delimiter="\t"))
    with mapping_path.open(encoding="utf-8", newline="") as handle:
        mappings = list(csv.DictReader(handle, delimiter="\t"))
    query_mapping = [
        row
        for row in mappings
        if row["raw_protein_id"] == query_raw_id and row["assembly_accession"] == "GCF_000007345.1"
    ]
    if len(query_mapping) != 1:
        raise ValueError(f"Expected one query mapping, observed {len(query_mapping)}")
    query_species = query_mapping[0]["species_id"]
    query_normalized_id = query_mapping[0]["normalized_protein_id"]
    query_fasta = normalized_directory / f"{query_species}.faa"
    query_sequences = _read_fasta(query_fasta)
    working_directory.mkdir(parents=True, exist_ok=True)
    single_query = working_directory / "query.faa"
    single_query.write_text(
        f">{query_normalized_id}\n{query_sequences[query_normalized_id]}\n",
        encoding="utf-8",
        newline="\n",
    )
    query_db = working_directory / "query_species.dmnd"
    _run([str(diamond_path), "makedb", "--in", str(query_fasta), "--db", str(query_db)])
    search = policy["search"]
    common = [
        "--outfmt",
        "6",
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "qlen",
        "slen",
        "evalue",
        "bitscore",
        "--more-sensitive",
        "--max-target-seqs",
        str(search["maximum_target_sequences"]),
        "--max-hsps",
        "1",
        "--evalue",
        str(search["evalue_maximum"]),
        "--threads",
        "4",
    ]
    all_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    commands: list[list[str]] = []
    selected = [row for row in panel if row["selection_status"] == "selected"]
    for species in sorted(selected, key=lambda row: row["species_id"]):
        target_species = species["species_id"]
        if target_species == query_species:
            continue
        target_fasta = normalized_directory / f"{target_species}.faa"
        target_db = working_directory / f"{target_species}.dmnd"
        makedb = [str(diamond_path), "makedb", "--in", str(target_fasta), "--db", str(target_db)]
        _run(makedb)
        commands.append(makedb)
        forward_path = working_directory / f"{target_species}.forward.tsv"
        forward = [
            str(diamond_path),
            "blastp",
            "--query",
            str(single_query),
            "--db",
            str(target_db),
            "--out",
            str(forward_path),
            *common,
        ]
        _run(forward)
        commands.append(forward)
        hits = _parse_hits(forward_path)
        target_sequences = _read_fasta(target_fasta)
        reciprocal_hits: dict[str, list[Hit]] = {}
        for hit in hits:
            reciprocal_query = working_directory / f"{target_species}.{hit.subject_id}.faa"
            reciprocal_query.write_text(
                f">{hit.subject_id}\n{target_sequences[hit.subject_id]}\n",
                encoding="utf-8",
                newline="\n",
            )
            reciprocal_path = (
                working_directory / f"{target_species}.{hit.subject_id}.reciprocal.tsv"
            )
            reciprocal_command = [
                str(diamond_path),
                "blastp",
                "--query",
                str(reciprocal_query),
                "--db",
                str(query_db),
                "--out",
                str(reciprocal_path),
                *common,
            ]
            _run(reciprocal_command)
            commands.append(reciprocal_command)
            reciprocal_hits[hit.subject_id] = _parse_hits(reciprocal_path)
        decision, reason, reciprocal_best, reciprocal_tie, paralog_count = classify_species(
            hits,
            reciprocal_hits,
            query_normalized_id=query_normalized_id,
            policy=policy,
        )
        best_tie = False
        if hits:
            tie_fraction = float(search["best_hit_tie_relative_bit_score"])
            best_tie = (
                sum(
                    hit.bit_score >= hits[0].bit_score * (1.0 - tie_fraction)
                    for hit in hits
                    if _passes(hit, search)
                )
                > 1
            )
        for rank, hit in enumerate(hits, 1):
            all_rows.append(
                {
                    "query_species": query_species,
                    "target_species": target_species,
                    "query_protein_id": query_raw_id,
                    "target_protein_id": hit.subject_id.split("__", 1)[-1],
                    "percent_identity": f"{hit.percent_identity:.6f}",
                    "alignment_length": hit.alignment_length,
                    "query_length": hit.query_length,
                    "target_length": hit.subject_length,
                    "query_coverage": f"{hit.query_coverage:.6f}",
                    "target_coverage": f"{hit.subject_coverage:.6f}",
                    "evalue": f"{hit.evalue:.8g}",
                    "bit_score": f"{hit.bit_score:.6f}",
                    "rank": rank,
                    "best_hit_tie": str(best_tie).lower(),
                    "reciprocal_best_hit": str(reciprocal_best if rank == 1 else False).lower(),
                    "reciprocal_tie": str(reciprocal_tie if rank == 1 else False).lower(),
                    "length_ratio": f"{hit.length_ratio:.6f}",
                    "fragment_flag": str(
                        hit.length_ratio
                        <= float(policy["classification"]["fragment_length_ratio_maximum"])
                        or hit.query_coverage < float(search["minimum_query_coverage"])
                        or hit.subject_coverage < float(search["minimum_subject_coverage"])
                    ).lower(),
                    "paralog_count": paralog_count,
                    "decision": decision,
                    "decision_reason": reason,
                    "rule_version": policy["policy_version"],
                }
            )
        summary_rows.append(
            {
                "target_species": target_species,
                "assembly_accession": species["assembly_accession"],
                "hit_count": len(hits),
                "top_target_protein_id": hits[0].subject_id.split("__", 1)[-1] if hits else "",
                "decision": decision,
                "decision_reason": reason,
                "best_hit_tie": str(best_tie).lower(),
                "reciprocal_best_hit": str(reciprocal_best).lower(),
                "reciprocal_tie": str(reciprocal_tie).lower(),
                "paralog_count": paralog_count,
            }
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(all_rows)
    summary_columns = (
        tuple(summary_rows[0])
        if summary_rows
        else (
            "target_species",
            "assembly_accession",
            "hit_count",
            "top_target_protein_id",
            "decision",
            "decision_reason",
            "best_hit_tie",
            "reciprocal_best_hit",
            "reciprocal_tie",
            "paralog_count",
        )
    )
    with summary_path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=summary_columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(summary_rows)
    decision_counts: dict[str, int] = {}
    for row in summary_rows:
        decision_counts[row["decision"]] = decision_counts.get(row["decision"], 0) + 1
    metadata = {
        "schema_version": "1.0",
        "tool": "DIAMOND",
        "tool_path": str(diamond_path.resolve()),
        "policy_path": str(policy_path.resolve()),
        "policy_sha256": _sha256(policy_path),
        "rule_version": policy["policy_version"],
        "query_species": query_species,
        "query_protein_id": query_raw_id,
        "comparison_species_count": len(summary_rows),
        "decision_counts": dict(sorted(decision_counts.items())),
        "commands": commands,
        "output_sha256": _sha256(output_path),
        "summary_sha256": _sha256(summary_path),
    }
    metadata_path.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return metadata


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diamond", required=True, type=Path)
    parser.add_argument("--normalized-directory", required=True, type=Path)
    parser.add_argument("--panel", required=True, type=Path)
    parser.add_argument("--mapping", required=True, type=Path)
    parser.add_argument("--policy", required=True, type=Path)
    parser.add_argument("--working-directory", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--metadata-output", required=True, type=Path)
    args = parser.parse_args()
    metadata = run_pilot(
        diamond_path=args.diamond,
        normalized_directory=args.normalized_directory,
        panel_path=args.panel,
        mapping_path=args.mapping,
        policy_path=args.policy,
        working_directory=args.working_directory,
        output_path=args.output,
        summary_path=args.summary,
        metadata_path=args.metadata_output,
    )
    print(json.dumps(metadata, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
