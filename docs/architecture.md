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
network, Excel, or external-service code. MVP-1A application services build the identifier
index, resolve queries, enumerate candidates, and apply configuration policies without running
biological evidence analysis. Ports describe capabilities. Local adapters parse input files,
while output adapters serialize already-validated models.

## Boundaries

- External services are isolated behind ports and have no implementation in MVP-0.
- JSONL is the canonical machine-readable evidence bundle.
- TSV is used only for flat records such as the manual structure queue.
- Excel is a derived schema view for human inspection.
- One parser or optional input failure is reported explicitly; missing evidence is not
  converted into a zero score.
- Local-only operation is the default and no network dependency is installed.

## Unimplemented engines

Gene context analysis, operon inference, orthology, phylogenetic profiles, domain analysis,
functional complementarity, localization prediction, known interaction retrieval, scoring,
tiers, and ranking remain outside MVP-1A. The pipeline emits `not_run` for each of these stages
and leaves every score and tier unset.

Identifier integration preserves original aliases, uses normalized values only for lookup,
and carries ambiguity forward instead of choosing a record. FASTA protein IDs are canonical
candidate identities.
