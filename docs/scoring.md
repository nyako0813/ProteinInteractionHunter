# Integrated scoring (MVP-1K)

Rule version: mvp1k-integrated-scoring-v1.

Scoring is disabled by default. It changes neither candidate disposition nor predicted relationship or evidence tier. Raw evidence records are not modified.

## Score formula

Only AVAILABLE components with a finite normalized value and nonzero configured weight are applied. Missing, not-run, not-applicable, failed, malformed, and zero-weight components remain in the breakdown but are excluded from the denominator.

Within each category, available configured weights are scaled proportionally when their sum exceeds the configured category cap. For each applied component:

    weighted_contribution = normalized_value * effective_weight
    raw_weighted_sum = sum(weighted_contribution)
    available_weight = sum(effective_weight)
    normalized_score = clamp(raw_weighted_sum / available_weight, 0, 1)
    provisional_score = normalized_score * output_scale

The formal output_score is emitted only when available_weight is at least minimum_evidence_weight and the number of represented categories is at least minimum_evidence_categories. Otherwise the provisional score is retained and rank is unset.

On the default 0–100 scale, 0 means no net support, 50 means half of the available capped weight supplies net support, and 100 means every applied component supplies maximum support. Absence from a database is neutral, not negative.

## Component normalization

| Component | Category | Positive | Neutral | Explicit negative | Unknown/excluded |
|---|---|---|---|---|---|
| genome_context | genomic_context | neighborhood with no intervening gene 0.8; neighborhood 0.5; other same-contig 0.1 | overlap or same feature 0 | none | missing coordinates, different-contig not-applicable |
| operon_proxy | genomic_context | supported 1.0; partial 0.3 | not supported 0 | none | unknown, missing, not-applicable |
| domain_pair | functional_annotation | any matched pair 1.0 | evaluated no-match 0 | reserved for an explicit biological-invalidity term | missing annotation |
| functional_complementarity | functional_annotation | matched rule 1.0 | evaluated no-match 0 | reserved for explicit functional contradiction | missing annotation |
| localization | cellular_compatibility | compatible 1.0 | evaluated but compatibility unknown 0 | incompatible negative contradictory-evidence penalty | missing annotation |
| orthology | evolutionary | supported 1.0, reduced by ambiguity penalty | evaluated unsupported 0 | none | missing orthology |
| phylogenetic_profile | evolutionary | similarity multiplied by shared-presence fraction | evaluated zero support 0 | at least 80% discordance with at least 3 informative species | insufficient informative species or missing profile |
| fusion | direct_interaction | supported 1.0 | no record or evaluated unsupported 0 | excessive component overlap | missing coverage |
| known_interactions | direct_interaction | direct 1.0; physical 0.9; functional 0.5; pair-only 0.4; predicted-only 0.2 | no record 0 | none in MVP-1K | missing or uncertain mapping |

Explicit negative strength is limited by penalties.contradictory_evidence. Candidate identifier ambiguity proportionally reduces all positive component values by penalties.ambiguous_mapping. Every application is recorded in applied_penalties.

## Categories and caps

The default groups are genomic context, functional annotation, cellular compatibility, evolutionary evidence, and direct interaction. Correlated components share a cap; available component weights are proportionally scaled, so adding a second correlated component cannot exceed that category cap.

## Ranking

Ranking is independent for each query. Only non-excluded candidates with sufficient evidence and a formal output score are eligible. Ordering is score descending, category count descending, available weight descending, then candidate ID ascending. Scores equal after tie_precision rounding receive the same dense rank (1, 1, 2). Secondary keys stabilize display order and do not break ties. Candidate TSV base row order is unchanged.

All arithmetic is performed in fixed component order and component, category, and aggregate values are rounded to 12 decimal places. Rank comparison separately uses the configured tie precision. Output and component breakdown preserve those deterministic values.
