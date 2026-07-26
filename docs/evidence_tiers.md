# Evidence tiers

Evidence tiers are a deterministic prioritization label for the strength, diversity,
specificity, and completeness of evidence available inside ProteinInteractionHunter. They are
not experimental confirmation, validation, or a claim that a pair truly interacts.

## Names and interpretation

- `tier_1`: strong, diverse evidence with at least one high-specificity association source.
- `tier_2`: comparatively strong support from multiple evidence categories.
- `tier_3`: limited but interpretable multi-category support.
- `tier_4`: evaluable evidence with weak or predominantly neutral support.
- `unclassified`: the pair is ineligible or does not meet the minimum quantitative gate.

Legacy A–X enum values remain readable for older artifacts, but the MVP-1L engine emits only the
five names above.

## Base tier

The engine uses the formal, already-rounded MVP-1K `output_score`; it never recalculates score or
rank. It tests tiers in the fixed order 1, 2, 3, 4 and assigns the first tier for which every
configured quantitative requirement passes:

| Tier | Score | Categories | Components | Available weight | High-specificity | Maximum negatives |
|---|---:|---:|---:|---:|---:|---:|
| 1 | >= 75 | >= 4 | >= 5 | >= 3.0 | >= 1 and required | 0 |
| 2 | >= 55 | >= 3 | >= 4 | >= 2.0 | >= 0 | 1 |
| 3 | >= 25 | >= 2 | >= 2 | >= 1.0 | >= 0 | 2 |
| 4 | >= 0 | >= 2 | >= 2 | >= 1.0 | >= 0 | 99 |

All defaults are configurable subject to monotonic validation: a lower tier cannot be stricter
than a higher tier. Thresholds must fit the configured scoring output scale.

## High-specificity evidence

Definition version `mvp1l-high-specificity-v1` recognizes these independent semantic classes:

- direct interaction support from known-interaction raw evidence;
- physical association support from known-interaction raw evidence;
- qualifying gene-fusion association evidence.

Multiple records from the same publication or source do not multiply these classes. Gene fusion
is an association signal, not proof of direct physical interaction. Functional association,
predicted interaction, co-expression, genome context, operon proxy, domain complementarity,
orthology, phylogenetic profile, and localization do not satisfy the Tier 1 high-specificity
gate by themselves.

## Caps

Caps are applied after the base tier in the fixed order below. A cap can only retain or lower a
tier; it never raises one.

1. `explicit_conflict`: an MVP-1K component with `applied=true` and `direction=negative`;
2. `predicted_only`: known-interaction types are predicted-only and no direct, physical, or
   fusion high-specificity signal exists;
3. `functional_association_only`: functional-association support exists without direct,
   physical, or fusion support.

The default cap for each condition is `tier_3`. Missing, `not_run`, `failed`,
`not_applicable`, neutral components, and an absent record are not conflicts.

## Eligibility and Unclassified

Formal classification requires the tier and scoring engines to be enabled, an available
integrated score, `sufficient_evidence=true`, a non-null formal score and rank, a non-self pair,
and a candidate that is not excluded. An enabled evaluation that fails any gate emits an
auditable `EvidenceTierResult` with `assigned_tier=unclassified` and the reason. When the engine
is disabled, the result list is empty, engine status is `not_run`, and the candidate tier remains
null.

Tier assignment does not change candidate inclusion, disposition, predicted relationship,
candidate order, score, rank, raw evidence, scoring components, or category breakdown. Tier is
not a ranking tie-breaker and never removes a candidate.

## Reproducibility

Rule version `mvp1l-evidence-tiers-v1`, thresholds, caps, the complete tier configuration,
high-specificity definition version, and the MVP-1K scoring rule dependency are recorded in the
run manifest. Base tier, assigned tier, semantic flags, applied caps, requirement outcomes,
support/conflict terms, warnings, and provenance are included in the canonical evidence bundle
and derived outputs.
