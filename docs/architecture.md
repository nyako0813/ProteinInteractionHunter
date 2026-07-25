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

## MVP-1B pipeline

```text
config + FASTA + GFF3 + optional annotation
  → identifier index and query resolution
  → all query-candidate pairs
  → normalized representative-feature index per contig
  → gene-context evidence per pair
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

## Gene-context boundary

MVP-1B reports coordinate facts only. It does not infer an operon, conserved neighborhood,
ortholog, functional relationship, physical interaction, rank, score, or Evidence Tier.
`gene_context` is available, missing, failed, not applicable, or not run independently for each
pair. All other evidence engines remain `not_run`; every score and tier remains `None`.

Linear contigs only are supported. Circular origin wrapping is not inferred. Contig lengths
come only from `##sequence-region`; when absent, edge distances are null and completeness is
`unknown`. Neighborhood membership uses representative feature-index distance and the configured
`gene_context.neighborhood_gene_count`, never a hidden base-pair threshold.

## Boundaries

- External services are isolated behind ports and have no MVP-1B implementation.
- JSONL is the canonical machine-readable evidence bundle.
- TSV and Excel are derived views for inspection.
- Parser or coordinate ambiguity is explicit; unavailable evidence is never converted to zero.
- Local-only operation is the default and no network dependency is installed.
- ProteinHunter_v5 is neither imported nor modified.

## Unimplemented engines

Operon inference, orthology, phylogenetic profiles, conserved-neighborhood analysis, domain
analysis, functional complementarity, localization prediction, fusion analysis, known
interaction retrieval, scoring, tiers, and ranking remain outside MVP-1B.
