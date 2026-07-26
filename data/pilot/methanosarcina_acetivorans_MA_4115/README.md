# MA_4115 minimal real-data pilot

## Scope

- Organism: *Methanosarcina acetivorans* C2A
- Assembly: GCF_000007345.1 (ASM734v1)
- Query protein accession: WP_011024006.1
- Current locus tag: MA_RS21465
- Old locus tags: MA_4115, MA4115
- GeneID: 1476009
- NCBI product annotation: alpha hydrolase
- Protein length: 200 aa
- Reference sequence: NC_003552.1, 5034570–5035172, minus strand

This minimal pilot checks candidate enumeration and local gene-context/operon-proxy evidence.
It is not an integrated interaction ranking and does not establish that MA_4115 experimentally
interacts with any candidate.

## Input provenance

The proteome and GFF are unmodified copies from the local NCBI dataset for
GCF_000007345.1 under ProteinHunter_v5. The query protein and nucleotide FASTA are unmodified
copies from the local `data/MA_4115_datasets/` NCBI query dataset.

| Pilot input | Source | SHA256 |
|---|---|---|
| `input/proteome.faa` | `ProteinHunter_v5/.../GCF_000007345.1/protein.faa` | `2eb0af39dece3e0c30bc673d4e10fab5888de221f8b4db6c71b658053062e6a2` |
| `input/genome.gff` | `ProteinHunter_v5/.../GCF_000007345.1/genomic.gff` | `d51451e44b3af446bb3f05a11c76986ca0895ca31f29a6ef92ab33b432c41368` |
| `input/query_MA_4115.faa` | `data/MA_4115_datasets/ncbi_dataset/data/protein.faa` | `17d137d509bf36d6379f917639319dde001c78a5c66772e5843f4984f640958a` |
| `input/query_MA_4115_gene.fna` | `data/MA_4115_datasets/ncbi_dataset/data/gene.fna` | `15a0c66704450e8ce15b815ebedd0dd7e0b62d78f32052726c7d6769faffd653` |

## Configuration

Only candidate generation, gene context, and its coupled operon proxy are enabled. Annotation,
domain, localization, orthology, phylogenetic-profile, fusion, and known-interaction tables are
not supplied. Scoring and evidence tiers are disabled because only the genomic-context evidence
category is available, while formal MVP-1 scoring requires at least two evidence categories.
Thresholds and weights are not relaxed for this pilot.

Paths in `config/pilot_minimal.yaml` are resolved relative to that YAML file.

## Run

```bash
cd /home/nyako/projects/ProteinInteractionHunter
.venv/bin/python -m protein_interaction_hunter validate-config \
  --config data/pilot/methanosarcina_acetivorans_MA_4115/config/pilot_minimal.yaml
.venv/bin/python -m protein_interaction_hunter validate-inputs \
  --config data/pilot/methanosarcina_acetivorans_MA_4115/config/pilot_minimal.yaml
/usr/bin/time -v .venv/bin/python -m protein_interaction_hunter generate-candidates \
  --config data/pilot/methanosarcina_acetivorans_MA_4115/config/pilot_minimal.yaml
```

- Pilot date: 2026-07-26
- Baseline commit: `8f9b41a6ed1ff2705cfa76dd858825c1aab866ce`

## Limitations

- No annotation table or external/local evidence table is used.
- FASTA descriptions and GFF identifiers are retained only where supported by current models.
- No annotation is inferred or supplemented.
- Scoring, ranking, evidence tiers, structure submission, and network access are disabled.
- Gene-context proximity and operon proxy support do not prove physical interaction.

## Pilot result

The minimal pilot completed successfully:

- 4,627 proteins and 4,627 query-candidate pairs
- 3,795 included, 831 fragment-policy flagged, and one excluded self pair
- 4,583 available gene-context records and 44 failed records caused by ambiguous multi-locus
  protein accessions
- one supported operon proxy: WP_011024007.1 (MA_RS21470), same minus strand, 5 bp separation,
  and no intervening gene
- scoring, rank, and evidence tiers remained not run
- wall time 12.34 seconds; peak resident memory 273,580 KiB

Audit helper outputs are `output/minimal/MA_4115_nearest_candidates.tsv` and
`output/minimal/MA_4115_genomic_neighborhood.tsv`. These are not canonical pipeline outputs.

## Annotation-table follow-up

The current annotation loader requires these exact TSV columns:

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

Additional provenance fields such as `old_locus_tag`, numeric `gene_id`, `description`,
coordinates, source version, assembly accession, and organism are useful for future preparation,
but are not accepted as substitutes for the required columns without a model change.

## NCBI annotation pilot

`scripts/build_ncbi_annotation_table.py` deterministically builds
`input/annotation_ncbi.tsv` from the matching NCBI RefSeq proteome FASTA and GFF3. The formal
loader table uses exactly these columns:

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

GFF CDS product is preferred; a FASTA description with only its final organism suffix removed
is the fallback when GFF product is absent. Gene name is populated only from an explicit `gene`
attribute. Functional category, localization, transmembrane topology, and confidence remain
empty because the NCBI inputs do not explicitly provide those fields in the required semantics.

The build found 4,596 formal rows for 4,627 proteins (99.330019%). All 4,627 FASTA accessions
occur in the GFF. Thirty-one multi-locus accessions are excluded from the formal one-row-per-ID
table rather than assigning an arbitrary locus. Same-locus split CDS records are merged.
`output/annotation/annotation_mapping_audit.tsv` retains every audited locus, and
`output/annotation/annotation_coverage.tsv` records coverage.

The annotation pilot differs from the minimal pilot only by run name, annotation input, and
output/cache/log directories. Gene context remains enabled; functional complementarity,
localization, scoring, and evidence tiers remain disabled. NCBI product text is not repurposed
as functional category or localization evidence. Scoring remains disabled because gene context
and operon proxy belong to the same genomic-context category.

```bash
cd /home/nyako/projects/ProteinInteractionHunter
.venv/bin/python scripts/build_ncbi_annotation_table.py \
  --fasta data/pilot/methanosarcina_acetivorans_MA_4115/input/proteome.faa \
  --gff data/pilot/methanosarcina_acetivorans_MA_4115/input/genome.gff \
  --annotation-output data/pilot/methanosarcina_acetivorans_MA_4115/input/annotation_ncbi.tsv \
  --audit-output data/pilot/methanosarcina_acetivorans_MA_4115/output/annotation/annotation_mapping_audit.tsv \
  --coverage-output data/pilot/methanosarcina_acetivorans_MA_4115/output/annotation/annotation_coverage.tsv \
  --annotation-source "NCBI RefSeq GCF_000007345.1" \
  --query-id WP_011024006.1
.venv/bin/python -m protein_interaction_hunter validate-config \
  --config data/pilot/methanosarcina_acetivorans_MA_4115/config/pilot_annotation.yaml
.venv/bin/python -m protein_interaction_hunter validate-inputs \
  --config data/pilot/methanosarcina_acetivorans_MA_4115/config/pilot_annotation.yaml
/usr/bin/time -v .venv/bin/python -m protein_interaction_hunter generate-candidates \
  --config data/pilot/methanosarcina_acetivorans_MA_4115/config/pilot_annotation.yaml
```

Next-stage evidence should come from explicit domain/function and localization/TM tools or
curated databases; no such external tool is run by this annotation pilot.

## InterProScan domain annotation follow-up

### Fixed software and environment

This workflow used official classic InterProScan 5 release `5.78-109.0` (released
2026-06-11; InterPro 109.0 and Pfam 38.2), downloaded from official EMBL-EBI FTP:

```text
https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.78-109.0/
interproscan-5.78-109.0-64-bit.tar.gz
MD5 3bb9a0794e9d69a0418a5298cdb04445
```

The verified archive was 7,037,746,512 bytes; the extracted installation is about 35 GiB at
`/home/nyako/tools/interproscan/5.78-109.0/` on the WSL Linux filesystem, not `/mnt/c`.
The runtime is Eclipse Temurin JRE `11.0.32+9` at
`/home/nyako/tools/java/temurin-11.0.32+9/`; its archive SHA256 is
`87ab4bf8dec10775d986957bc313816678f9227f1d033d7d6e6a1d00dace5b95`. No system
Python or repository environment was changed. The launcher heap maximum was reduced from 15 GiB
to 8 GiB because WSL has 15 GiB RAM and 4 GiB swap. The host exposed 16 logical CPUs and about
948 GiB was free on the Linux filesystem before installation.

Official documentation used:

- <https://interproscan-docs.readthedocs.io/en/v5/ReleaseNotes.html>
- <https://interproscan-docs.readthedocs.io/en/v5/InstallationRequirements.html>
- <https://interproscan-docs.readthedocs.io/en/v5/HowToRun.html>
- <https://interproscan-docs.readthedocs.io/en/v5/OutputFormats.html>
- <https://interproscan-docs.readthedocs.io/en/v5/ImprovingPerformance.html>

InterProScan software is Apache-licensed, while bundled member databases retain their own terms.
Only installed Pfam 38.2 was used. Licensed components such as SignalP, Phobius, and TMHMM were
not acquired or used. `-dp` disabled EBI precalculated-match lookup, so no sequences were
uploaded. Current InterPro mappings plus GO and pathway fields were requested with
`-goterms -pa`.

### Commands and raw runs

Installed `interproscan.sh --help` was checked first. The compact equivalents of the timed
commands are:

```bash
interproscan.sh -i query_MA_4115.faa -f TSV -appl Pfam -cpu 2 -dp   -goterms -pa -vtsv -b MA_4115_pfam
interproscan.sh -i interproscan_subset.faa -f TSV -appl Pfam -cpu 4 -dp   -goterms -pa -vtsv -b subset_pfam
interproscan.sh -i proteome.faa -f TSV -appl Pfam -cpu 4 -dp   -goterms -pa -vtsv -b full_pfam
```

All exited 0. Query-only took 55.43 s (peak RSS 3,764,908 KiB). The deterministic subset had 55
proteins/31,430 aa and took 66.28 s (peak RSS 3,619,488 KiB). Selection comprised the query,
unique ±10 CDS, minimal-pilot nearest 20, the first eight accession-sorted hypothetical or
uncharacterized products, the first eight products matching membrane/transporter/permease/
symporter/antiporter, and the eight longest and eight shortest proteins. IDs were deduplicated
and final FASTA accession-sorted. Membership and reasons are in
`output/interproscan/subset_selection_audit.tsv`.

The subset produced 166 Pfam hits across 40/55 proteins (72.727273%); the other 15 were valid
no-hits, not errors. It had 53 unique Pfam accessions, 39 proteins with InterPro accessions, 21
with GO terms, 25 with pathways, 12 multi-hit proteins, and 8 repeated-signature proteins. Raw
TSV size was 77,792 bytes.

The environment was therefore `READY_TO_RUN` for conservative Pfam-only full proteome.
The 4,627-protein run completed in 5:41.86 (peak RSS 3,689,820 KiB), exit 0, no swaps, and a
6,727,739-byte raw TSV. UTC times, commands, version, stdout/stderr, and `time -v` are under
`output/interproscan/`. This conclusion does not promise the resources required by all
bundled analyses.

### Converter and formal table

`scripts/convert_interproscan_domains.py` validates InterProScan 5 TSV, accepts comments and
empty files, preserves IDs/analysis/signature/description/coordinates, sorts deterministically,
and writes UTF-8/LF. Only exact duplicate input rows are removed. Non-overlapping repeats remain
separate; a conflicting formal identity is fatal. Bad column counts or coordinates are fatal.
Proteome-unknown IDs remain in audit but are excluded from formal output. The audit preserves
InterPro accession/description, score or e-value, status/date, GO, pathways, MD5/length, version,
and source path.

The formal loader requires exactly:

```text
protein_id  source  accession  name  start  end  architecture_index
```

`name` may be empty per row; identity fields and valid 1-based coordinates are required.
There are no formal score, e-value, status, InterPro, GO, or pathway columns. The loader rejects
malformed/duplicate identities, retains overlaps and repeats, and accepts a header-only table.
A domain table alone is not interaction evidence: a separately sourced domain-pair rule YAML and
domains on query and candidate are also required.

The full conversion produced 5,618 formal rows for 3,421/4,627 proteins (73.935595%), 1,206
no-hits, and 1,909 unique Pfam accessions. InterPro coverage was 3,336 proteins (72.098552%), GO
1,668 (36.049276%), and pathways 2,002 (43.267776%). It found 1,138 multi-hit and 301 repeated-
signature proteins. Of 1,085 hypothetical/uncharacterized proteins, 49 had a hit; 3,372/3,542
other proteins had a hit. Exact duplicates, unknown IDs, malformed rows, and excluded rows were
all zero.

### MA_4115 result and interpretation

| Protein | Analysis | Signature | Description | Start-end | Score/e-value | InterPro | GO | Pathway |
|---|---|---|---|---:|---|---|---|---|
| WP_011024006.1 | Pfam | PF24167 | Family of unknown function (DUF7411) | 1-194 | 2.7E-76 | IPR055834, Protein of unknown function DUF7411 | none | none |

The single signature covers 194/200 aa (0.97), with no overlap or repeat. It supports DUF7411
membership but does not identify a catalytic family. It neither independently validates nor
directly contradicts broad NCBI `alpha hydrolase` annotation because Pfam/InterPro labels
remain function-unknown. It cannot establish substrate, mechanism, tRNA binding, sulfur
insertion, cnm5U involvement, or interaction partner.

### Coverage-only pilot and readiness

`config/pilot_domains.yaml` points at the converted table and retains gene context plus NCBI
annotation. `domains.enabled` is intentionally false and `rules_path` null because no
real, versioned, biologically sourced MA_4115 domain-pair rules exist. Synthetic fixture rules
were not reused and domain names were not turned into rules manually. Validation loaded 5,618
rows on 3,421 known proteins, with zero unknown IDs and one query hit.

The timed run produced the same 4,627 pairs, order, disposition, flags, annotation, gene context,
and operon proxy as annotation-only. Excluding `run_id`, candidate TSV SHA256 was identical:
`84c2c606cc8768271de057021d3e988fc231d3e1df02117572f976de4956340d`. Warning TSV was
also identical: `3ba5b25a021000b79e24fe151f4e8c51beef2ae493160658bc7941b469a922b3`.
The manifest records the domain table and SHA256. Domain evidence, scores, ranks, and tiers remain
empty/not run. Runtime was 11.85 s and peak RSS 281,796 KiB.

Domain-pair readiness is `NOT_READY`: annotation coverage is adequate, but there is no
curated rule source/version, known-pair derivation, or biological basis for a formal match.
Functional readiness is `NOT_READY`: GO/InterPro coverage exists, but no reviewed mapping to
the functional engine vocabulary and no real complementarity rules exist. Mapping/provenance
policy is a separate future stage; this workflow does not auto-fill `functional_category`.

Raw InterProScan outputs, full converted table, subset FASTA, logs, cache, and large audits stay
ignored. Intended commit candidates are converter, tests, config, README, and minimal validation/
manifest changes. Work started at `f2b2c1208d7765328c1a0ee8c9144c9cf3200246`; no commit or
push was performed.
