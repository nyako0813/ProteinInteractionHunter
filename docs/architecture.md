# Architecture

## Dependency direction

```text
Domain models
    ↑
Application services
    ↑
Ports
    ↑
Adapters and outputs
```

Domain models contain scientific records and validation rules only. They do not import CLI,
network, Excel, or external-service code. Application services build the identifier index,
resolve queries, enumerate candidates, normalize GFF coordinates, and calculate observation-only
gene context. Ports describe capabilities. Local adapters parse input files, while output
adapters serialize already-validated models.

## MVP-1 pipeline

```text
config + local FASTA/GFF/annotation/evidence tables
  → identifier index and all query-candidate pairs
  → independent raw evidence engines
  → category-capped integrated scoring and per-query ranking (optional)
  → evidence tier classification (optional; requires scoring)
  → canonical JSONL + derived TSV/Excel + manifest/warnings
```

The GFF loader retains 1-based closed coordinates, source, feature type, strand (`+`, `-`, or
`?`), identifiers, all parents, and `##sequence-region` boundaries. The coordinate index uses
FASTA protein IDs as canonical identities. For a protein, one unique CDS coordinate is the
representative; a related gene may have a wider extent and is retained as a warning. Multiple
distinct CDS coordinates are ambiguous and are never selected by row order. Coordinates on
contradictory replicons are fatal.

Gene/CDS rows belonging to one protein or gene unit collapse to one representative feature.
Unmapped gene/RNA units remain observable representatives, preventing gene/CDS double counting
while preserving non-protein features for `intervening_feature_count`. Features are sorted by
`(start, end, representative_id)` for deterministic indices.

## Evidence boundaries

Each engine reports available, missing, failed, not applicable, or not run independently.
Disabled or unavailable evidence is excluded from the scoring denominator and is never converted
to negative evidence. Operon is an explicit proxy, gene fusion is association evidence rather
than proof of direct physical interaction, and Evidence Tier is a prioritization label rather
than experimental confirmation.

Linear contigs only are supported. Circular origin wrapping is not inferred. Contig lengths
come only from `##sequence-region`; when absent, edge distances are null and completeness is
`unknown`. Neighborhood membership uses representative feature-index distance and the configured
`gene_context.neighborhood_gene_count`, never a hidden base-pair threshold.

## Boundaries

- External services remain isolated behind ports and are not enabled by MVP-1.
- JSONL is the canonical machine-readable evidence bundle.
- TSV and Excel are derived views for inspection.
- Parser or coordinate ambiguity is explicit; unavailable evidence is never converted to zero.
- Local-only operation is the default and no network dependency is installed.
- ProteinHunter_v5 is neither imported nor modified.

## Current engine scope

All configured MVP-1 local evidence, scoring, ranking, and tier engines are implemented.
`_UNIMPLEMENTED_ENGINES` is empty. Conserved-neighborhood analysis, network retrieval, and
automatic structure prediction remain outside the current local pipeline.
