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
network, Excel, or external-service code. Application services coordinate validation without
performing candidate analysis. Ports describe capabilities. Local adapters parse fixture
files, while output adapters serialize already-validated models.

## Boundaries

- External services are isolated behind ports and have no implementation in MVP-0.
- JSONL is the canonical machine-readable evidence bundle.
- TSV is used only for flat records such as the manual structure queue.
- Excel is a derived schema view for human inspection.
- One parser or optional input failure is reported explicitly; missing evidence is not
  converted into a zero score.
- Local-only operation is the default and no network dependency is installed.

## Unimplemented engines

Candidate generation, gene context analysis, operon inference, orthology, phylogenetic
profiles, domain analysis, functional complementarity, localization prediction, known
interaction retrieval, scoring, tiers, and ranking remain outside MVP-0. The pipeline
placeholder raises an explicit exception rather than returning empty or synthetic results.
