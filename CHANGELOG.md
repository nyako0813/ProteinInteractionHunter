# Changelog

## Unreleased - MVP-1A

- Added versioned identifier normalization and ambiguity-preserving alias resolution.
- Added multiple-query resolution, deterministic duplicate groups, fragment/hypothetical
  classification, and explicit candidate disposition policies.
- Added initial all-`not_run` evidence bundles with uncomputed scores.
- Added `generate-candidates`, candidate TSV, effective config snapshot, warning summary, and
  completed run provenance.

## 0.1.0 - 2026-07-25

- Initialized an independent Python 3.12 `src`-layout project.
- Added strict configuration, domain models, ports, local fixture loaders, and output writers.
- Added versioned JSON Schemas with drift validation.
- Added synthetic FASTA, GFF3, annotation, and expected-output fixtures.
- Added CLI validation commands and an explicitly unimplemented analysis pipeline.
- Added unit, integration, separation, lint, type-check, and schema-validation foundations.
