# MA_4115 formal scoring integration

## Decision

The MA_4115 real-data pilot is ready for biological review under the existing
`mvp1k-integrated-scoring-v1` and `mvp1l-evidence-tiers-v1` contracts. The score, dense
rank, and tier are prioritization outputs from currently available computational evidence.
They are not evidence that a physical interaction occurs.

Readiness at the audited start commit `5662d3deb4936f1a0c38593c41b9743c76d2a3a7`:

| Item | Decision |
|---|---|
| Gene context formal evidence | READY |
| Operon proxy formal evidence | READY |
| Localization formal evidence | READY |
| Orthology formal evidence | READY |
| Phylogenetic profile formal evidence | READY |
| Functional complementarity | NOT_READY |
| Domain-pair evidence | NOT_READY |
| Integrated scoring | READY |
| Rank | READY |
| Evidence tier | READY |
| Final candidate list | READY_FOR_BIOLOGICAL_REVIEW |

Functional complementarity remains disabled: the query category is
`trna_modification`, but there are zero query-applicable complementarity rules. Domain
coverage is annotation-only and is not domain-pair evidence.

## Inputs and validation

All files were regular, readable files. The FASTA contained 4,627 unique proteins; the
query `WP_011024006.1` occurred once. There were no unknown IDs in the domain, orthology,
or profile tables, no duplicate identities, and no malformed rows. Forty-four proteins
lacked unambiguous coordinates and 31 lacked annotation rows.

| Input | Rows or records | Size (bytes) | SHA-256 |
|---|---:|---:|---|
| `input/proteome.faa` | 4,627 proteins | 1,752,636 | `2eb0af39dece3e0c30bc673d4e10fab5888de221f8b4db6c71b658053062e6a2` |
| `input/genome.gff` | 9,973 lines; 4,583 mapped proteins | 2,911,965 | `d51451e44b3af446bb3f05a11c76986ca0895ca31f29a6ef92ab33b432c41368` |
| `input/annotation_ncbi.tsv` | 4,596 data rows | 420,873 | `f8ccdc5e50a430592823331a5e93df40e289986f7b3261deabcdab10cfad4d2b` |
| `input/domains_interproscan.tsv` | 5,618 domains; 3,421 proteins | 384,161 | `e034c4d0a3b3bcc3309cb06531ba975bda5f80cddee5124f9e0d7bb04de52453` |
| `input/annotation_psortb.tsv` | 4,596 data rows | 578,560 | `e202a41b7ed273efeac634d4a07cd6ad46efb0bda1d908b9cf2a1182d638c1c8` |
| `input/orthology_orthofinder.tsv` | 54,467 records; 4,056 proteins | 9,945,406 | `c47f0cea2ac43d643e42700b9f1dda240ca6b0f7214bc84a74c7d5b63ac6d4bc` |
| `input/phylogenetic_profiles.tsv` | 120,302 cells; 4,627 × 26 | 14,711,889 | `7447e4668f1ad4e1efa85a43ef1947edcdff23518df4255644290d021d6a6295` |
| `input/annotation_functional.tsv` | 4,596 data rows | 578,577 | `46b4ee5e53c5546f9fa17a01fd243ff4a98915b80fcc5e71b87380f3cb5cc3e9` |

The query values were: 200 aa, locus `MA_RS21465`, product `alpha hydrolase`,
localization `cytosolic`, category `trna_modification`, domain `PF24167`, 22 orthology
rows, and 26 profile rows. Table schemas matched the repository converters. Relevant
source/rule provenance was NCBI RefSeq `GCF_000007345.1`, PSORTb
`3.0.6+dfsg-3build4` with `psortb-3.0-archaea-terse-mapping-v1`, OrthoFinder `3.1.5`
with `ma4115-orthofinder-orthology-v1`, profile build
`ma4115-orthology-profile-v1`, formal profile rule `mvp1h-phylogenetic-profile-v1`,
and the versioned functional-category mapping. The domain and functional tables were
validated and included in provenance, but their formal engines stayed disabled.

## Enabled and disabled engines

Enabled:

- gene context and its coupled operon proxy;
- annotation-only localization compatibility;
- local-table orthology;
- local-table phylogenetic profile;
- integrated scoring and evidence tiers.

Disabled:

- functional complementarity, domain pair, known interactions, fusion, structure
  evidence, DeepTMHMM, SignalP, and automatic structure prediction.

No hidden score component was observed. Every disabled component serialized as
`not_run`, had effective weight zero, and was excluded from the denominator.

## Scoring contract

Components group into exactly three available categories in this pilot:

| Category | Components |
|---|---|
| `genomic_context` | `genome_context`, `operon_proxy` |
| `cellular_compatibility` | `localization` |
| `evolutionary` | `orthology`, `phylogenetic_profile` |

Configured weights are 1.0, 1.0, 0.5, 0.75, and 1.0 respectively. Disabled component
weights remain in the configuration but do not enter the available denominator. Category
caps are 1.5, 0.5, and 2.0 for the three categories above (with unchanged 1.5 and 2.0
caps for the disabled functional and direct categories). The contradiction and ambiguity
strengths are 0.25 and 0.10.

Missing/not-run/not-applicable evidence is excluded from the denominator. Available
neutral evidence contributes its effective weight and zero numerator. Explicit
localization incompatibility and strong profile discordance contribute `-0.25` before
weighting. The final numerator is divided by available weight, clamped to `[0, 1]`, and
scaled to `[0, 100]`. A formal output needs weight at least 1.0 and at least two
categories. Ties use score rounded to eight decimals. Dense rank is assigned by rounded
score, with category count, available weight, and stable candidate ID as deterministic
sort keys. The dense rank changes only when the rounded score changes.

An excluded self pair retains an auditable internal score trace but receives no rank and
an `unclassified/self_pair` tier. Flagged non-self candidates remain rank eligible under
the existing contract.

## Single-engine pilots

Every pilot retained 4,627 candidates and the baseline dispositions: 3,795 included, 831
flagged, and one excluded self. Scoring and tiers remained ungenerated.

| Pilot | Result |
|---|---|
| Gene context | 4,583 available; 44 failed coordinate mappings |
| Operon proxy | 1 supported; 2,282 partial; 2,299 unsupported; 44 unknown; 1 not applicable |
| Localization | 3,838 available; 789 missing; 2,891 compatible; 947 incompatible |
| Orthology | 4,056 available; 571 missing; 1 supported; 4,055 unsupported |
| Profile | 4,627 available; 586 supported; 4,041 unsupported |

The sole operon-supported pair was `WP_011024007.1`: 5 bp, same minus strand, and no
intervening gene. This is a proxy support signal, not a known positive interaction.
Gene proximity likewise does not establish a physical interaction.

Localization mapped 2,891 candidates to cytosolic, 840 to membrane, and 107 to secreted;
789 were missing. Unknown or unsupported localization values remained neutral/missing,
not incompatible. PSORTb confidence was not converted to an additional score.

Profile similarity ranged from 0 to 1 for 4,625 candidates; two profiles with no
informative value were excluded from the component denominator. The exact 0.8 threshold
was inclusive (13 candidates). The nearest observed values were 0.7916666667 below and
0.8076923077 above. Ambiguous and fragment-only observations remained unknown, not
absence/presence.

## Category-cap audits

For 4,582 non-self coordinate-resolved pairs, the two configured genomic weights summed
to 2.0 and were scaled to 0.75 each, exactly reaching the 1.5 category cap. The self pair
had one applicable genomic component (weight 1.0); 44 coordinate-failed pairs had none.
Thus two genomic components count as one category.

The evolutionary configured total was at most 1.75, below its 2.0 cap, so the cap was
non-binding in this dataset. This is expected at the unchanged formal weights. Grouping
still prevents orthology and profile from inflating the category count. Support
combinations were: profile only 585, both 1 (the excluded self), orthology only 0, and
neither 4,041.

## Integrated evidence and score

Available category counts were:

| Categories | Candidates |
|---:|---:|
| 3 | 3,825 |
| 2 | 771 |
| 1 | 31 |
| 0 | 0 |

Available weights ranged from 1.0 to 3.75. Formal scores were generated for 4,596
candidates, including the excluded self's trace; 31 had `None` from the minimum-category
rule. Of the generated scores, 4,247 were positive and 349 were zero.

| Statistic | Value |
|---|---:|
| Ranked non-self candidates | 4,595 |
| Minimum | 0 |
| Minimum positive | 0.0701754386 |
| Maximum | 71.8974358974 |
| Mean | 19.2724731507 |
| Median | 18.3333333333 |
| Q25 / Q75 | 9.1111111111 / 28.5128205128 |
| Q90 / Q95 / Q99 | 37.7435897436 / 39.7948717949 / 43.8974358974 |
| Unique rounded scores | 677 |
| Tie groups | 456 |
| Largest tie | 349 |

Negative component counts were 0 for 2,685 candidates, 1 for 1,747, and 2 for 195.
There were 2,137 contradiction penalty records and no ambiguity penalty records.

The independent audit recomputed raw component strength from evidence records, configured
weight, category assignment, cap scaling, weighted contribution, available weight,
numerator, normalized/final score, eight-decimal tie value, and dense rank. All 41,643
component comparisons across 4,627 candidates matched; total failures were zero.

## Rank determinism

The scoring-only run was repeated unchanged and repeated with deterministic reversal of
the annotation/localization, orthology, and profile table row orders. Candidate order,
score, and rank were identical. Full JSONL was identical after removing only `run_id`.
TSV, JSONL, and Excel score/rank values agreed in all three runs. Run-specific IDs,
paths, and timestamps were not compared as scientific output.

## Evidence tiers

The unchanged tier thresholds produced:

| Tier | Candidates |
|---|---:|
| Tier 1 | 0 |
| Tier 2 | 3 |
| Tier 3 | 1,543 |
| Tier 4 | 3,049 |
| Unclassified | 32 |

Unclassified comprises 31 insufficient-category pairs and the excluded self. The maximum
score (71.8974) is below the Tier 1 score threshold and the pilot has only three available
categories and no high-specificity direct/physical/fusion evidence.

The current `mvp1l-high-specificity-v1` contract defines `predicted_only` narrowly:
known-interaction types must be non-empty and contain only `predicted`, with no
high-specificity signal. Because known interactions are disabled, `predicted_only=false`
for all 4,627 pairs. Consequently, the predicted-only Tier 3 cap is correctly
non-operative under the current definition, and three quantitatively eligible pairs
retain Tier 2. These are `WP_011024007.1` (71.8974), `WP_011024005.1` (57.8974), and
`WP_011024009.1` (56.0). This is an important interpretation constraint: Tier 2 here
does not mean experimentally supported evidence.

Explicit conflict was present for 1,942 candidates. All already had base Tier 3/4 or
were unclassified, so the Tier 3 explicit-conflict cap did not change an assigned tier.
The cap condition and order were nevertheless serialized and independently inspected.
No cap was silently removed.

Observed score boundaries were: no score at exactly 75 or 55; nearest values below were
71.8974 and 51.8974; nearest above 55 was 56.0. There was no exact 25; nearest values were
24.9797570850 and 25.0303030303. Tier 4's exact zero boundary contained 349 candidates.
The machine-readable tier decision audit preserves every satisfied/failed requirement and
applied cap for all candidates.

## Top candidates and biological sanity review

`WP_011024007.1` ranked first with score 71.8974358974 and Tier 2. It had gene-context
strength 0.8, operon support 1.0, cytosolic compatibility, neutral unsupported orthology,
and profile component strength 0.8461538462. It is not treated as a known correct pair.

The next unique dense ranks included:

| Rank | Candidate | Product | Score | Tier |
|---:|---|---|---:|---|
| 1 | `WP_011024007.1` | DNA-binding protein | 71.8974 | 2 |
| 2 | `WP_011024005.1` | 50S ribosomal protein L39e | 57.8974 | 2 |
| 3 | `WP_011024009.1` | adenylosuccinate synthase | 56.0000 | 2 |
| 4 | `WP_011024001.1` and four tied neighbors | local ribosomal/translation proteins | 51.8974 | 3 |
| 5 | `WP_248698120.1` | 50S ribosomal protein L18Ae | 50.8718 | 3 |
| 6 | three tied candidates | thermosome/PGK/TBP | 48.0000 | 3 |

Dense rank 25 had score 41.8461538462 and a 30-member tie, so the Top 25 review contains
185 candidates. Dense rank 50 had score 38.4444444445 and one candidate, giving 297
candidates in the Top 50 review. Boundary ties are included, never truncated.

In the Top 50 review set (including ties), all 297 were cytosolic and had orthology
annotations, 247 passed the profile support thresholds, 26 were fragment-flagged, one was
hypothetical, and none was membrane/secreted. Review-only keyword counts were: ribosomal
17, translation factors 7, metabolic-enzyme patterns 192, transcription regulators 7,
transporters 1, RNA-related 33, sulfur-related 4, ATPase 5, and redox-related 18. These
keywords were never passed to scoring.

No self, excluded candidate, unknown protein ID, or score-`None` candidate was ranked.
Repeated dense ranks correspond only to exact eight-decimal score ties; gaps follow the
dense-rank contract.

## Baseline A/B and warnings

Against `output/functional_coverage`, candidate ID sets and all selected invariant TSV
fields were identical. Candidate JSON, raw gene-context evidence, operon evidence, and
disposition were also identical. There were zero unexpected TSV or JSON differences.
Expected changes were activation of localization/orthology/profile statuses, scoring,
rank, and tier outputs plus their sheets and manifest fields.

All baseline warnings were retained: 44 coordinate-mapping warnings (also summarized as
44 gene-context failures), 831 fragment flags, 1,085 hypothetical/uncharacterized flags,
and 31 missing annotations. New expected summaries were 31 insufficient-score warnings,
789 missing localization annotations, 571 missing orthology annotations, and profile
unknown-observation summaries. Component traces may carry the same coordinate warning on
both gene-context and operon components, but the candidate/run warning summary deduplicates
the cause to one count per candidate.

## Performance and output

All runs used one worker, exited zero, swapped zero pages, and left raw artifacts ignored.

| Run | Wall time | Peak RSS |
|---|---:|---:|
| Gene context | 14.55 s | 285,180 KiB |
| Localization | 4.29 s | 232,956 KiB |
| Orthology | 23.58 s | 908,380 KiB |
| Profile | 9.30 s | 410,132 KiB |
| Orthology + profile | 26.04 s | 1,117,912 KiB |
| Integrated evidence, no score | 42.05 s | 1,261,288 KiB |
| Scoring only | 50.83 s | 1,610,764 KiB |
| Scoring repeat | 55.97 s | 1,618,240 KiB |
| Scoring with reversed rows | 52.19 s | 1,610,152 KiB |
| Tier-enabled final | 58.94 s | 1,693,560 KiB |

Final output sizes were 17,495,831 bytes TSV, 140,965,242 bytes JSONL, 11,932,598 bytes
Excel, 3,154 bytes warning summary, 11,984 bytes manifest, and 4,716 bytes config
snapshot. None is tracked by Git.

## Reproduction and artifacts

The tracked formal configuration is
`data/pilot/methanosarcina_acetivorans_MA_4115/config/pilot_formal_scoring.yaml`.
Run it with:

```bash
.venv/bin/python -m protein_interaction_hunter generate-candidates \
  --config data/pilot/methanosarcina_acetivorans_MA_4115/config/pilot_formal_scoring.yaml
```

Use `scripts/audit_ma4115_formal_scoring.py` against the baseline, single-engine,
scoring-repeat, shuffled, and final output directories. It emits ignored
`summary.json`, `recalculation_failures.json`, `tier_decision_audit.tsv`, and Top 25/50
review tables.

## Limitations and next precision stage

The top list is dominated by computational co-localization, neighborhood, and
co-evolution signals. Shared absence is part of the existing profile metric, orthology
and profile originate from the same OrthoFinder panel, localization is PSORTb-only, and
fragment/hypothetical flags remain common. None provides direct physical-interaction
evidence.

The next precision stage should add independently curated, versioned, query-applicable
functional/domain-pair rules; validated direct or physical interaction records; expanded
orthology/profile panels with ambiguity audits; and, only after independent validation,
structure or experimental evidence. Thresholds, weights, caps, and tier rules must not be
fit retrospectively to the current top candidates.
