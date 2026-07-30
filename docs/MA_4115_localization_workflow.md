# MA_4115 localization and transmembrane-evidence workflow

## Scope and decision

This workflow prepares auditable real-data localization input for
*Methanosarcina acetivorans* C2A query WP_011024006.1 (MA_RS21465; old locus
MA_4115/MA4115). It does not infer localization from product, Pfam, or InterPro
text. PSORTb, DeepTMHMM, and SignalP have distinct roles:

- PSORTb: predicted subcellular-localization source.
- DeepTMHMM: predicted alpha-helical transmembrane-topology source.
- SignalP: predicted N-terminal signal-peptide source.

These sources are not interchangeable observations. A SignalP positive is not
automatically extracellular; one TM helix is not automatically a membrane
protein; PSORTb Unknown is not cytosolic; no signal peptide is not cytosolic;
and no TM helix is not automatic soluble-interaction support. No multi-source
consensus rule or duplicate scoring is introduced here.

The implementation is deliberately limited to a PSORTb converter, deterministic
small-pilot subset builder, full PSORTb run, and coverage-only pipeline audit.
Localization evidence, integrated scoring, and evidence tiers remain disabled.

## Existing formal contract audit

`LocalizationConfig` has only `enabled: bool` and
`source: annotation_only`; there is no independent localization-table adapter.
The formal route is `LocalAnnotationTsvLoader`, with exactly these columns:

```text
protein_id
gene_name
locus_tag
product
functional_category
localization_annotation
transmembrane_annotation
annotation_source
annotation_confidence
```

The loader requires one row per protein ID. It cannot retain multiple
localization sources independently, and its one confidence and source field
apply to the entire annotation row. The generic evidence model can retain
multiple provenance records with method and metadata, but the annotation loader
cannot populate that richer source-specific representation.

The current annotation-only engine normalizes free text and recognizes
`s_layer`, `secreted`, `membrane`, and `cytosolic`. Candidate/query compartment
equality is compatible, inequality is contradictory, and an unrecognized or
missing compartment has compatibility `None`. Explicit transmembrane free text
can represent `none`, `single_pass`, `multi_pass`, or generic
`transmembrane`; the model has an integer helix count but the formal loader has
no coordinate or orientation columns. Signal-peptide text can become true,
false, or missing, but signal type, cleavage site, and probabilities have no
formal columns.

`EvidenceStatus` is `available`, `missing`, `not_applicable`, `failed`, or
`not_run`; origin can distinguish annotation from local prediction. If scoring
were enabled, localization is the independent `cellular_compatibility`
category: compatible is positive, explicit incompatibility receives the
configured contradiction penalty, compatibility `None` is neutral, and
missing/not-run is excluded. It is not `genomic_context`.

Validation loads the annotation table even when localization is disabled. The
manifest fingerprints it as an input. JSONL stores localization evidence in
the candidate bundle; candidate TSV exposes localization status, compartments,
compatibility, topology, source, and confidence; Excel has a
`Localization_Evidence` sheet. With localization disabled these evidence
records are absent/not-run.

Because merely adding the 31 proteins absent from the existing 4,596-row NCBI
annotation table would change candidate annotation-presence state in a
coverage-only run, the converter preserves the base table's ID scope. All
4,627 PSORTb predictions remain in raw/audit coverage, while predictions for
the 31 base-table-absent IDs are audit-only. This avoids a core schema change
and keeps the A/B candidate contract stable.

PSORTb `Cellwall` is also audit-only: it is not mapped to `s_layer`, because
generic archaeal cell-wall prediction does not establish an S-layer and
`cell_wall` is not in the current formal compatibility vocabulary. Adding a
reviewed formal cell-wall compartment is a separate schema/vocabulary decision.

## External-tool readiness

### PSORTb

Official PSORTb documentation identifies version 3.0, Archaea mode `-a`, whole
FASTA input, and normal/terse/long outputs. Terse output is header plus
`SeqID`, `Localization`, and `Score`; `Unknown` means no site passes the 7.5
cutoff. The official download page lists Bio-Tools-PSort 3.0.6, recommends the
Docker distribution until a newer release, supports Ubuntu/Linux, and states
GNU GPL. Sources:

- <https://psort.org/documentation/index.html>
- <https://psort.org/downloads/>
- <https://github.com/brinkmanlab/psortb-docker>

The installed runtime is isolated at:

```text
/home/nyako/tools/psortb/3.0.6-ubuntu24.04/
```

Ubuntu 24.04 package `psortb 3.0.6+dfsg-3build4` and dependencies were
downloaded and extracted there without system installation. The PSORTb package
SHA256 is
`b9ee84206e41515cca3f5c49b5bc5bff4d1aab5b76df097e5b6a7b48fea319b2`;
its package copyright records GPL-2+ for the main code, GPL-2+ for ModHMM,
and Artistic for the local SVM wrapper. Extracted runtime size is about
238 MiB. The system Python and repository venv were not changed.

Installed `psort --help` was the option source of truth and installed
`psort --version` reported `PSORTb version 3.0`. Ubuntu's
`ncbi-blast+-legacy 2.12.0+ds-4build2` wrapper was added in the same isolated
root (package SHA256
`80708d5c25ba4a5ada853a871a4c27a56b0b8623f44fa678c5ce931fb40b5568`).
The packaged archaeal SCL-BLAST FASTA was indexed locally with BLAST+
`makeblastdb`; this reported a few invalid residues and one unusually long
title in the upstream database, and therefore remains a reproducibility
warning. A query-only dependency check then exited 0 without PSORTb/BLAST
warnings.

Resource expectation measured on WSL Ubuntu (16 logical CPUs, 15 GiB RAM) is
single-process CPU execution. Query-only took 0.68 s and 85,308 KiB peak RSS;
the 82-protein subset took 39.94 s and 118,164 KiB peak RSS.

Status: `READY`.

### DeepTMHMM

The official DTU DeepTMHMM 1.0 service predicts alpha-helical and beta-barrel
classes plus residue-level topology, covers all domains of life, and is stated
to scale to proteomes. DTU directs local download through BioLib. The hosted
software is available for academic and commercial use, but commercial
on-premises use requires a commercial license. The public Docker image is
large (about 4.4 GiB compressed), old, tagged only `latest`, and does not give
this pilot an approved model/version pin. No licensed/approved local
distribution is present, so no unofficial package is substituted and no
sequence is uploaded.

Source: <https://services.healthtech.dtu.dk/services/DeepTMHMM-1.0/>

Status: `BLOCKED_PENDING_LICENSED_SOFTWARE`.

When approved software is available, a separate converter should retain
sequence-level TM status, helix count, and topology class formally only if the
loader contract is reviewed. Residue segments, start/end, inside/outside,
signal-like labels, model version, and raw topology remain audit fields.

### SignalP

The official SignalP 6.0 service supports Archaea and predicts Sec/SPI,
Sec/SPII, Tat/SPI, Tat/SPII, Sec/SPIII, Other, cleavage sites, region labels,
and probabilities. Fast mode is recommended for batches; slow mode is about
six times slower and is intended for accurate region boundaries. Server
downloads include JSON, one-line prediction summary, processed FASTA, and
GFF3. The portable Python package is offered through an academic download
route; other users must contact the DTU software manager. No approved local
package/license is present, so the service was not used and no proteome was
uploaded.

Sources:

- <https://services.healthtech.dtu.dk/services/SignalP-6.0/>
- <https://services.healthtech.dtu.dk/services/SignalP-6.0/7-Portable.php>

Status: `BLOCKED_PENDING_LICENSED_SOFTWARE`.

SignalP is only a signal-peptide source. A future audit must retain protein ID,
class/type, cleavage site, probability, mode, and version without directly
mapping a positive result to a final localization.

## Input verification

The source proteome and GFF are readable regular files and match the expected
checksums:

```text
protein.faa  2eb0af39dece3e0c30bc673d4e10fab5888de221f8b4db6c71b658053062e6a2
genomic.gff  d51451e44b3af446bb3f05a11c76986ca0895ca31f29a6ef92ab33b432c41368
```

The proteome contains 4,627 unique proteins and exactly one
WP_011024006.1 query. Existing ignored pilot copies are byte-identical, so no
tracked duplicate is added. Since PSORTb returns the FASTA title as SeqID, the
runtime FASTAs use accession-only headers while preserving sequences:

```text
full PSORTb FASTA SHA256:
38285255eeb636d1a7cc4a16085ebd9bf27489d6038502470cbb8f475194ceb1

subset PSORTb FASTA SHA256:
2c85e5ec2a4440c7f24f0626494b0ddba78c18f54405c47c220aee3e4b7557b0
```

## Converter contract

`scripts/convert_psortb_localization.py` supports official PSORTb 3 terse
output only. This avoids guessing the version-specific long format. It:

- accepts the exact optional header, comments, blank lines, and header-only
  input;
- requires exactly three columns and a finite score from 0 through 10;
- preserves SeqID verbatim;
- rejects conflicting duplicate IDs and audits exact duplicate rows;
- audits unknown and missing proteome IDs and requires the requested query;
- sorts formal and audit outputs deterministically and writes UTF-8/LF;
- records source/version/role/command/path, Archaea mode, raw row number,
  mapping decision, exclusions, and raw/formal checksums;
- never infers from product, Pfam, InterPro, signal peptide, or TM text.

Mapping rule `psortb-3.0-archaea-terse-mapping-v1` is explicit:

| Raw PSORTb | Formal localization | Status |
|---|---|---|
| Cytoplasmic | cytosolic | accepted |
| CytoplasmicMembrane | membrane | accepted |
| Cellwall | empty | unsupported/audit-only |
| Extracellular | secreted | accepted |
| Unknown | empty | unknown/missing |

The formal output has exactly the existing nine annotation columns. Base NCBI
rows and product/locus fields are preserved. PSORTb source is appended at row
level, but its 0–10 score is not placed in the formal 0–1 annotation-confidence
field; score and exact source metadata remain in audit. `transmembrane_annotation`
is never populated by PSORTb.

The audit table has:

```text
protein_id raw_prediction raw_score formal_localization mapping_status
mapping_rule_version known_protein_id included_in_formal_table duplicate_kind
raw_row_number normalization_decision exclusion_reason source_name
source_version source_role organism_mode source_command source_file
raw_input_sha256 formal_output_sha256
```

Coverage distinguishes raw rows, known proteins represented, non-Unknown,
Unknown, absent, malformed, duplicate, unknown-ID, unsupported, excluded, and
localization distribution. Metadata separately fingerprints formal, audit, and
coverage outputs.

## Small pilot

`scripts/build_psortb_subset.py` starts from the prior deterministic real-data
subset and adds accession-sorted representatives for multiple Pfam domains, no
Pfam hit, non-standard residues if present, highest N-terminal hydrophobic
fraction, and lowest maximum 19-residue hydrophobic fraction. It emits
accession-only FASTA and a reasoned audit. The final 82 unique proteins include:

- query, 20 unique ±10 CDS/neighborhood entries, and operon-supported
  WP_011024007.1;
- 20 deterministic nearest candidates;
- eight membrane/transporter/permease products;
- eight multiple-Pfam and eight no-Pfam representatives;
- eight longest and eight shortest valid proteins;
- eight N-terminal-hydrophobic and eight no-obvious-hydrophobic representatives.

No non-standard residues occur in this proteome. IDs are deduplicated and
accession-sorted.

Exact compact PSORTb command:

```bash
psort -a -o terse \
  -r /home/nyako/tools/psortb/3.0.6-ubuntu24.04/root/usr/lib/psort \
  -m /home/nyako/tools/psortb/3.0.6-ubuntu24.04/root/usr/bin \
  psortb_subset.faa
```

The run exited 0 in 39.94 s, peak RSS 118,164 KiB, with no PSORTb/BLAST
warnings. It returned 82 rows: Cytoplasmic 47, CytoplasmicMembrane 15,
Cellwall 2, Extracellular 1, and Unknown 17. Duplicate, malformed, unknown-ID,
and missing subset-ID counts were zero. Raw output SHA256 was
`b4bd7dabaac1ddfa614bbf88f51ff88caa879c792584a21a71ff77e310adbd48`.
The converter and formal loader accepted it.

MA_4115 was returned verbatim as:

```text
WP_011024006.1    Cytoplasmic    7.50
```

It maps to formal `cytosolic`; no TM topology or signal peptide was inferred.

## Full-proteome result

Exact compact command:

```bash
psort -a -o terse \
  -r /home/nyako/tools/psortb/3.0.6-ubuntu24.04/root/usr/lib/psort \
  -m /home/nyako/tools/psortb/3.0.6-ubuntu24.04/root/usr/bin \
  psortb_proteome.faa
```

The full local run completed successfully:

- exit status: `0`; output: 4,628 lines (header plus 4,627 proteins)
- elapsed wall time: `34:06.78`; maximum RSS: `182,280 KiB`
- raw output: 156,903 bytes; SHA-256
  `c0b6fdfe0f6c2436698c68b5e90e17be887642640b4eacc242aab0f7c38d5c12`
- raw distribution: Cytoplasmic 2,920; CytoplasmicMembrane 840;
  Extracellular 107; Cellwall 30; Unknown 730 (15.776961%)
- converter audit: 4,627 unique known proteins, with zero missing proteins,
  unknown IDs, malformed rows, exact duplicates, or conflicting duplicates
- formal table: 4,596 rows, 3,838 accepted mapped predictions; 730 Unknown,
  30 unsupported Cellwall, and 31 predictions outside the preserved base-table
  scope remain explicit in the audit rather than becoming negative evidence
- formal table SHA-256:
  `e202a41b7ed273efeac634d4a07cd6ad46efb0bda1d908b9cf2a1182d638c1c8`
- `WP_011024006.1`: Cytoplasmic, score 7.50, formally mapped to `cytosolic`

## Coverage-only pipeline audit

`config/pilot_localization_coverage.yaml` points at the enriched formal
annotation table. It preserves the prior domain coverage-only configuration:

```text
domains.enabled: false
localization.enabled: false
scoring.enabled: false
evidence_tiers.enabled: false
```

Thus validation and the manifest read/fingerprint the PSORTb-enriched table,
but no localization evidence reaches candidates or scoring.

Both `validate-config` and `validate-inputs` succeeded. Candidate generation
completed with exit status 0 in 13.87 seconds (maximum RSS 285,956 KiB) and
produced 4,627 query-candidate pairs: 3,795 included, 831 flagged, and one
excluded self-candidate.

A structured A/B comparison against the prior `output/domains` coverage-only
run passed. After removing only `run_id`, candidate TSV rows and JSONL evidence
bundles were identical in content and order; warning summaries were byte
identical. Every candidate had an empty localization list, localization,
scoring, and evidence-tier statuses of `not_run`, null score/rank/tier values,
and empty integrated-scoring/evidence-tier traces. The run manifest fingerprints
`annotation_psortb.tsv` with SHA-256
`e202a41b7ed273efeac634d4a07cd6ad46efb0bda1d908b9cf2a1182d638c1c8`.

## Readiness summary

```text
PSORTb installation: READY
PSORTb raw prediction: READY
PSORTb conversion: READY
PSORTb formal localization input: READY_FOR_COVERAGE_ONLY
DeepTMHMM: BLOCKED_PENDING_LICENSED_SOFTWARE
SignalP: BLOCKED_PENDING_LICENSED_SOFTWARE
Localization engine pilot: COVERAGE_ONLY_READY
Formal localization evidence: NOT_ENABLED_PENDING_REVIEW
Formal scoring: NOT_READY
```

Before formal localization evidence is enabled, the project must decide how
to represent cell wall, multiple sources, source-specific confidence, and
future TM/signal-peptide observations without double counting. PSORTb alone
can create false confidence, and same predicted compartment is compatibility
support rather than proof of interaction. Query or candidate Unknown must stay
missing/neutral. DeepTMHMM or SignalP and orthology/phylogenetic-profile
evidence are still absent. Therefore score, rank, and tier must remain unset.
