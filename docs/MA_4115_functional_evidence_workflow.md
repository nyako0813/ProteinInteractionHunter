# MA_4115 functional annotation and evidence-readiness workflow

## Scope and decision

This workflow extends the Methanosarcina acetivorans C2A MA_4115 real-data pilot to functional annotation. It does not claim a physical interaction. The run began from clean main at 7d6812af1c84c24ba4e12490010636a9f3ed2d17, equal to origin/main. The query is WP_011024006.1, current locus MA_RS21465, old locus MA_4115.

| Layer | Status | Reason |
|---|---|---|
| Raw eggNOG annotation | READY | Exact local run, schema and ID audits passed |
| Formal category mapping | READY for query; sparse proteome-wide | One exact curated KO mapping |
| MA_4115 query category | READY | ko:K07585 maps to trna_modification |
| Versioned complementarity rules | NOT_READY | No supported query/candidate role rule set |
| Functional coverage-only input | READY | Loader and pipeline validation passed |
| Formal functional evidence | NOT_READY | Zero query-applicable rules and evaluable pairs |
| Formal scoring | NOT_READY | Disabled; no score, rank, or tier was emitted |

No description keyword, NCBI product, Pfam or InterPro label, COG category, shared pathway, or EC number was converted directly to interaction evidence.

## Query and input audit

Source proteome:

    /home/nyako/projects/ProteinHunter_v5/data/databases/target/methanosarcina_acetivorans/ncbi_dataset/data/GCF_000007345.1/protein.faa

It is a readable regular file with SHA256 2eb0af39dece3e0c30bc673d4e10fab5888de221f8b4db6c71b658053062e6a2, 4,627 records, no duplicate IDs, and no empty sequences. WP_011024006.1 occurs exactly once, is 200 aa, and has sequence SHA256 1a9d9b85ef933cfc0aea94bf8ac95bc9f51894f97b4ec5262bc1d0fa49b7ddd2. Its original header is WP_011024006.1 alpha hydrolase [Methanosarcina acetivorans]. The exact extracted query FASTA has SHA256 17d137d509bf36d6379f917639319dde001c78a5c66772e5843f4984f640958a.

## Existing engine contract

LocalAnnotationTsvLoader requires exactly these columns:

    protein_id
    gene_name
    locus_tag
    product
    functional_category
    localization_annotation
    transmembrane_annotation
    annotation_source
    annotation_confidence

There is one row per protein and duplicate IDs are fatal. There is no separate functional-annotation loader. functional_category is one optional string; the core loader defines no multi-category delimiter, vocabulary registry, or unknown-category validation. The pipeline combines product, functional category, and FASTA description into normalized annotation text. Source terms therefore remain in a separate audit. The builder preserves all NCBI and localization columns and only replaces functional_category.

The current engine is keyword-role based, not a formal category-pair engine. Its rule model has ruleset version, role IDs with include/exclude terms, and directional query-role to candidate-role pair rules. Symmetry requires an explicit reverse rule. Missing query or candidate text yields MISSING. An evaluated non-match is AVAILABLE with matched false and no_matching_pair_rule, not negative evidence. Multiple matching roles and pair rules can yield multiple records.

The schema has no reference, taxonomic scope, support level, explicit positive/negative/neutral relation, wildcard, priority, or duplicate/conflict validation. conflicting_terms are labels, not a formal contradictory rule type. The schema is sufficient for a disabled coverage-only connection, but not for production MA_4115 rules. No core schema was changed and synthetic fixture rules were not reused.

When enabled, functional complementarity contributes to scoring component functional_complementarity, grouped under category functional_annotation with unchanged cap 1.5. Scoring still requires at least two evidence categories and uses the existing missing-value denominator policy. Functional, scoring, and tiers can remain disabled with rules_path null. This reuses the optional-input audit pattern from the other pilots.

## Source readiness

### eggNOG-mapper

[eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) stable v2 was selected. Exact source is [tag v2.1.15](https://github.com/eggnogdb/eggnog-mapper/tree/v2.1.15), tag object 3e35ae46a37f175b7555ab0eaa88aac841b09452, peeled commit 74cec65609153afd0ffea41957f988c263907a59, dated 2026-05-22. Software license is AGPL-3.0. v3 remained testing or unreleased with an incompatible DB and was not used.

Installation is isolated at /home/nyako/tools/eggnog-mapper/2.1.15; source is under source, venv under .venv, Python is 3.12.3, and total size is 349 MiB. Upstream Biopython 1.76 cannot build on Python 3.12, so the controlled compatibility environment uses Biopython 1.87 and NumPy 2.5.1 while retaining psutil 5.7.0 and XlsxWriter 1.4.3. Installed source is editable. The executable reports emapper-v2.1.15, expected DB 5.0.2, and bundled DIAMOND 2.0.11.

Reproducible installation equivalents are git clone at tag v2.1.15, python3 -m venv for the isolated environment, pip install of Biopython 1.87, NumPy 2.5.1, psutil 5.7.0 and XlsxWriter 1.4.3, followed by pip install -e of the fixed source checkout. The fixed tag and commit above are the source of truth.

The DB came from the official [emapperdb-5.0.2 directory](https://eggnog5.embl.de/download/emapperdb-5.0.2/) because the old downloader hostname no longer resolved. All runs were local DIAMOND searches with no upload. Local archive birth times were 2026-07-30 13:42:13 UTC for eggnog.db.gz, 14:02:04 UTC for taxa, and 14:02:37 UTC for DIAMOND. GO files were created locally at 14:23:39 to 14:23:40 UTC. The 59 GiB DB remains outside Git. eggNOG data is CC BY-NC 4.0 and is not redistributed.

| File | Bytes | SHA256 |
|---|---:|---|
| eggnog.db.gz | 6,776,977,123 | 20fd3e57d96ab19f2094131e859cb0f2e634e10f594c01ef34335076317a8e0a |
| eggnog.taxa.tar.gz | 72,797,584 | ead5f9f1cb42a51a990c3980fcb272a2f6c07b41c75654739e799d96a190a7f3 |
| eggnog_proteins.dmnd.gz | 5,208,806,170 | c74fe142c6d0f1c749a086f754f1658874c0a52d0927626b9d263098c136b320 |
| eggnog.db | 41,370,988,544 | 4d6857c129423889edcd607995b28d97944787a588224f31dd7f745fe9ef1e4f |
| eggnog.taxa.db | 278,003,712 | 84d56ba1cbb091cc4c4064ba0f8235020654ed50fd6d617f829c4e45551faffd |
| eggnog.taxa.db.traverse.pkl | 6,628,719 | de7bb2b35489c9e5e1f2522575c8107ebd3d7117f8910dc6745aab68ac2d2173 |
| eggnog_proteins.dmnd | 9,285,439,161 | f4b63db0313021e4a2ca1d47a6b3f7862a3c8c2016abf3b51ef5865e653642a7 |

### InterProScan and COG2024

Existing Pfam 38.2 and InterPro 109.0 results were reused without overwriting the formal domain table. They contain 5,618 Pfam hits on 3,421 proteins, InterPro assignments on 3,336, GO on 1,668, and pathways on 2,002. These transferred annotations are not independent evidence from eggNOG.

The COG2024 update was published on 2025-01-06 and its current incremental files were updated through August 2025. The current [NCBI COG service](https://www.ncbi.nlm.nih.gov/research/cog/) covers 2,103 bacterial and 193 archaeal species. The [COG2024 data index](https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2024/data/) defines protein-to-COG, definition, category, pathway, organism, and taxonomy files. Official MD5 values include cog-24.cog.csv 6e0ba7bbd09eb422a545f5c75fff4355, cog-24.def.tab 1df3d05b4ca26144bc9a975f21ccd695, cog-24.fun.tab 47e95258577aa11b5840f4431e8ce66e, cog-24.mapping.tab c1269bacaf9a180adfb1b9e0c805a912, cog-24.pathways.tab 2c066a903def5bd032db47f2fa71c165, COGorg24.faa.gz 664b003570a2dc0b0268870711cab17e, and COGorg24.gene.tab.gz 6d9cebd88f4c9a08fd7fb3a9ea6747aa.

The official API directly returns WP_011024006.1 from GCF_000007345.1, locus MA_RS21465, as a full-length membership-class-0 match to COG2117, bit score 281 and E-value 1.25e-97. COG2117 is Adenosine-binding enzyme, AANH-like superfamily, category J, with no COG pathway. This agrees with eggNOG COG assignment but is broader than the KO result. It was audited, not added as a second formal source.

### Gene Ontology

Fixed folder is /home/nyako/databases/go/2026-06-19; OBO data-version is releases/2026-06-15.

| File | Bytes | SHA256 |
|---|---:|---|
| go-basic.obo | 32,215,811 | c72fc198a86983d55e43aac585d1ffdbeb6e3601475b3f18b6045acdc0a0734c |
| go-computed-taxon-constraints.obo | 1,488,560 | c7c3422b00ae02fde76a0222af49dc24c386a58c87db9fcd708dfbb2d7766adc |

Sources are official [ontology downloads](https://geneontology.org/docs/download-ontology/), [release folder](https://release.geneontology.org/2026-06-19/ontology/), and [taxon constraints](https://geneontology.org/docs/taxon-constraints/). Only is_a, and policy-authorized part_of, may propagate. regulates is not inherited. MF and BP remain separate; CC is excluded. Roots, unknown IDs, and obsolete IDs are excluded. replaced_by and consider values remain in audit reasons. Every ancestor path is retained deterministically and cycles are fatal.

The eggNOG file does not expose per-assignment evidence codes or NOT qualifiers. The run used go_evidence non-electronic, and this limitation is preserved. No GO mapping was accepted in v1, so no taxon-constrained GO mapping could activate. The constraint file is pinned for a future mapping version; taxon-conflict enforcement remains an explicit readiness gap.

## Runs and commands

All commands used protein input, DIAMOND, local DB files, tax_scope archaea, tax_scope_mode inner_narrowest, target_orthologs all, go_evidence non-electronic, pfam_realign none, report_orthologs, and report_no_hits. E-value and seed thresholds were 0.001. Query and subset used DIAMOND ctg; full used auto. The full command is preserved in output converter metadata as:

    /home/nyako/tools/eggnog-mapper/2.1.15/.venv/bin/emapper.py -i data/pilot/methanosarcina_acetivorans_MA_4115/input/proteome.faa --itype proteins -m diamond --cpu 8 --data_dir /home/nyako/databases/eggnog/5.0.2 --dmnd_algo auto --sensmode sensitive --dmnd_iterate yes --evalue 0.001 --seed_ortholog_evalue 0.001 --tax_scope archaea --tax_scope_mode inner_narrowest --target_orthologs all --go_evidence non-electronic --pfam_realign none --report_orthologs --report_no_hits --output full --output_dir data/pilot/methanosarcina_acetivorans_MA_4115/output/functional/eggnog/full --temp_dir data/pilot/methanosarcina_acetivorans_MA_4115/output/functional/eggnog/full

| Run | Proteins or aa | Wall | Peak RSS KiB | Exit | Raw disk |
|---|---:|---:|---:|---:|---:|
| Query | 1 / 200 | 23.43 s | 2,323,780 | 0 | 52 KiB |
| Subset | 104 / 48,997 | 37.09 s | 2,635,740 | 0 | 252 KiB |
| Full | 4,627 | 4:37.00 | 3,514,028 | 0 | 8.1 MiB |

All runs had zero swaps. UTC times, stdout, stderr, exact commands, and time -v logs remain in ignored run directories.

Deterministic subset SHA256 is 25fc098ecc73e4aea6b03a48b8867d8763898855c182fb5646513bb8acaf12a3; membership audit SHA256 is 0f70c9cd4732d67b9e661d615bf2738903d3e830c200e2899492b575646a002b. Its fixed policy covers the query, plus or minus 10 CDS, WP_011024007.1, orthology/profile selections, paralog and unassigned cases, InterPro GO/Pfam presence and absence, localization strata, functional product strata, and length extremes. It is accession-sorted and deduplicated.

Subset results: 95 hits and seeds, 90 annotations, five hit-without-annotation records, nine no-hits, and zero missing, unknown, duplicate, or malformed IDs. Coverage: description 84, preferred name 31, COG 84, GO MF 0, BP 0, CC 2 with 12 terms total, EC 20, KO 52, pathway 32, module 23, reaction 14, PFAM 75.

## Raw converter and full coverage

scripts/convert_eggnog_annotations.py validates the exact 21-column v2.1 header, metadata comments, row widths, IDs, duplicates, query preservation, hit/seed/no-hit status, deterministic order, UTF-8/LF, versions, command, and checksums. It produces separate audit, coverage, and JSON metadata. Fields include protein ID, seed and search statistics, eggNOG OGs, COG, description, preferred name, GO, EC, KEGG fields, BRITE, TC, CAZy, BiGG, PFAMs, version, command, raw row, status, and exclusion reason.

| Full metric | Count | Percent |
|---|---:|---:|
| Proteome | 4,627 | 100 |
| Search hit or seed | 4,059 | 87.724227 |
| Annotated row | 3,642 | 78.711908 |
| Description | 3,221 | 69.613140 |
| Preferred name | 665 | 14.372163 |
| eggNOG OG | 3,642 | 78.711908 |
| COG | 3,221 | 69.613140 |
| GO MF proteins | 23 | 0.497082 |
| GO BP proteins | 22 | 0.475470 |
| GO CC proteins | 16 | 0.345796 |
| GO unknown-aspect proteins | 5 | 0.108061 |
| Total GO assignments | 911 | n/a |
| EC | 859 | 18.564945 |
| KEGG KO | 1,839 | 39.744975 |
| KEGG pathway | 1,074 | 23.211584 |
| KEGG module | 775 | 16.749514 |
| KEGG reaction | 678 | 14.653123 |
| PFAM field | 3,076 | 66.479360 |
| No hit | 568 | 12.275773 |
| Hit, no annotation | 417 | 9.012319 |
| Unknown-function description | 250 | 5.403069 |
| Broad-only description | 130 | 2.809596 |

Malformed, duplicate, unknown, and unexplained missing IDs were all zero. Raw checksums: annotations 3412ea01ef1963b68fa4b2316fee511924304f19c75639eededcea7cc477b300; hits cab611b413d8772e73dc9a038bd845ae8ff9974fa359392afc6a3edd17c9637a; seed a5638fa4ddb6a43b05db37428b8ebb4604e75cd7df4b1f4de4676efe85471b33; orthologs 30cb7727cdcc903463f4266be45b397d37f37815c59f3987e40a282e0a3dc069; full audit 62ca90b45893069c6ed6607e66a2b238f51bb9ad9d4196e44491937bb0584553. Full conversion rerun outputs were byte-identical.

## MA_4115 and source union

MA_4115 status is annotated. Seed is 192952.MM_0804, E-value 3.69e-128, score 365.0, 91 percent identity, and full query/target coverage. OGs are COG2117, arCOG00037, 2XTGQ, and 2N9S4; max level is Methanomicrobia; COG category is J. Description is subunit of tRNA(5-methylaminomethyl-2-thiouridylate) methyltransferase contains the PP-loop ATPase domain. KEGG KO is ko:K07585. Preferred name, GO, EC, pathway, module, reaction, and PFAM are empty.

Existing query InterPro result is PF24167/IPR055834 DUF7411 with no GO/pathway. eggNOG provides orthology-transferred description and KO; InterPro provides the domain. Neither overwrites the other. Proteome union: eggNOG GO 36 proteins, InterPro GO 1,668, exact shared GO 14, ancestor/descendant six, unresolved 34, conflict-like zero. eggNOG PFAM occurs on 3,076 and InterPro Pfam on 3,421; exact shared accession is zero because eggNOG uses names while InterPro formal output uses accessions. InterPro remains domain source of truth.

## Formal vocabulary and mapping

Versioned policy is data/pilot/methanosarcina_acetivorans_MA_4115/config/functional_category_mapping.v1.yaml, version ma4115-functional-category-mapping-v1. It defines only trna_modification. The category has explicit definition, criteria, sources, aspects, relations, taxonomic applicability, evidence quality, ambiguity, and conflict policies.

The only accepted rule is exact ko:K07585 to trna_modification, priority 100. References are curated [KEGG K07585](https://rest.kegg.jp/get/K07585), which defines a tRNA methyltransferase and lists MA_4115, and [mac:MA_4115](https://rest.kegg.jp/get/mac:MA_4115), which assigns the same 200-aa protein to K07585. This is not inferred from alpha hydrolase, DUF7411, eggNOG description, COG J, or ATPase text.

Policy validation rejects schema deviations, duplicates, unknown category, invalid status/match/relation, empty rationale/reference, GO roots, obsolete or unreachable GO, aspect mismatch, wrong version, priority ties, and conflicting accepted mappings.

Formal output preserves 4,596 existing annotation records; 31 multi-locus accessions intentionally lack a formal row. Counts: one formal protein, zero multiple categories, one exact mapping, zero ancestor, zero ambiguous, zero conflict, 3,783 only-unmapped/source-no-category, and 843 no-source. Query is the one formal protein. Annotation SHA256 is 46b4ee5e53c5546f9fa17a01fd243ff4a98915b80fcc5e71b87380f3cb5cc3e9; mapping audit is 29f2a2f88f19efd7e62fb813f55139c61d7efb42fb490c0d02d15fcc87874891. Independent builds were byte-identical.

## Complementarity readiness

No complementarity rule file was created. A curated query identity does not establish complementary candidate roles, physical partners, or contradictions. Shared KO, GO, COG, pathway, or EC is correlated annotation, not a justified rule. Orthology produced eggNOG transfer and is also existing evolutionary input, so these signals are not independent.

Readiness counts: one formal category; one query category; zero candidate categories beyond the query; zero rules; zero query-applicable rules; zero evaluable/supported/contradictory/neutral pairs; 4,626 non-self pairs unevaluated. False-positive risk is high if broad J, ATPase, or pathway similarity is promoted. False-negative risk remains from sparse mappings and absent pair biology. Formal functional evidence is therefore NOT_READY.

## Coverage-only pipeline and A/B

pilot_functional_coverage.yaml connects enhanced annotation and existing domain, localization, orthology, and profile inputs for validation/manifest audit. functional_complementarity is disabled with null rules; scoring and tiers are disabled.

Validation loaded 4,627 proteins, 4,583 coordinates, 5,618 domains, 54,467 orthology records, and 120,302 profile observations without unknown optional IDs. Run result: 4,627 pairs, 3,795 included, 831 flagged, one excluded self-pair. Runtime 11.00 s, peak RSS 286,000 KiB, swaps zero.

Against orthology/profile baseline, candidate TSV and JSONL excluding run ID, warning summary, candidate order/count, and Excel sheet names/dimensions were identical. Functional, score, rank, and tier remained not_run or empty. Expected manifest changes were command, timestamps, config path/hash/snapshot, recorded git commit, input metadata, run ID, and run name. Disposition, fragment flags, coordinates, gene context, operon, annotation, domain, localization, orthology, and profile evidence did not change.

## Limitations and next step

- eggNOG annotations are orthology transfers correlated with existing orthology evidence.
- GO transfer is sparse and lacks per-assignment evidence code and NOT qualifier detail.
- Taxon constraints are pinned but not enforced because v1 accepts no GO mapping. Future GO mapping must implement and test this first.
- Core functional rules lack production provenance, polarity, scope, priority, and conflict semantics.
- Annotation loader has no formal multi-category contract. v1 exercises only one category on one protein.
- No score, rank, tier, shadow score, threshold relaxation, or automatic evidence was produced.

The next stage is literature or curated reaction/complex review for a versioned role-pair rule candidate audit with explicit scope and polarity. If that cannot be done without speculation, formal functional evidence must remain NOT_READY and scoring must remain disabled.
