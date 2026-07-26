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
