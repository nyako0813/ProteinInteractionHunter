"""Strict YAML configuration models and config-relative path resolution."""

from pathlib import Path
from typing import Literal

import yaml
from pydantic import Field

from protein_interaction_hunter.exceptions import ConfigurationError
from protein_interaction_hunter.models.base import StrictModel


class ProjectConfig(StrictModel):
    name: str = Field(min_length=1)
    run_name: str = Field(min_length=1)
    random_seed: int = 0
    local_only: bool = True


class InputConfig(StrictModel):
    proteome_fasta: Path
    genome_gff: Path
    annotation_table: Path | None = None
    reference_proteomes: list[Path] = Field(default_factory=list)


class QueryConfig(StrictModel):
    protein_ids: list[str] = Field(min_length=1)
    allow_multiple: bool = True

    @classmethod
    def _non_empty_ids(cls, values: list[str]) -> list[str]:
        normalized = [value.strip() for value in values]
        if any(not value for value in normalized):
            raise ValueError("query protein IDs must not be empty")
        if len(set(normalized)) != len(normalized):
            raise ValueError("query protein IDs must be unique")
        return normalized

    def model_post_init(self, __context: object) -> None:
        self.protein_ids = self._non_empty_ids(self.protein_ids)
        if not self.allow_multiple and len(self.protein_ids) > 1:
            raise ValueError("multiple query IDs are disabled")


class CandidateGenerationConfig(StrictModel):
    include_hypothetical_proteins: bool = True
    include_missing_coordinates: bool = True
    minimum_length_aa: int = Field(default=20, ge=1)
    duplicate_sequence_policy: Literal["flag", "exclude", "error"] = "flag"
    fragment_policy: Literal["flag", "exclude", "include"] = "flag"
    self_candidate_policy: Literal["exclude"] = "exclude"


class GeneContextConfig(StrictModel):
    enabled: bool = False
    neighborhood_gene_count: int = Field(default=5, ge=0)
    require_query_coordinates: bool = False
    operon_proxy_max_intergenic_bp: int = Field(default=200, ge=0)


class OrthologyConfig(StrictModel):
    enabled: bool = False
    source: Literal["computed", "local_table"] = "computed"
    engine: Literal["blast", "diamond"] = "blast"
    database: Path | None = None
    local_table: Path | None = None


class PhylogeneticProfileConfig(StrictModel):
    enabled: bool = False
    source: Literal["local_table"] = "local_table"
    local_table: Path | None = None
    minimum_shared_species: int = Field(default=2, ge=1)
    minimum_informative_species: int = Field(default=3, ge=1)
    minimum_profile_similarity: float = Field(default=0.8, ge=0.0, le=1.0)


class FusionConfig(StrictModel):
    enabled: bool = False
    source: Literal["local_table"] = "local_table"
    local_table: Path | None = None
    minimum_supporting_records: int = Field(default=1, ge=1)
    minimum_component_coverage: float = Field(default=0.6, ge=0.0, le=1.0)
    maximum_component_overlap_fraction: float = Field(default=0.2, ge=0.0, le=1.0)


class EnabledConfig(StrictModel):
    enabled: bool = False


class DomainConfig(StrictModel):
    enabled: bool = False
    source: Literal["none", "local_table"] = "none"
    local_table: Path | None = None
    rules_path: Path | None = None


class FunctionalComplementarityConfig(StrictModel):
    enabled: bool = False
    rules_path: Path | None = None


class LocalizationConfig(StrictModel):
    enabled: bool = False
    source: Literal["annotation_only"] = "annotation_only"


InteractionType = Literal[
    "physical",
    "direct",
    "genetic",
    "functional_association",
    "co_complex",
    "co_expression",
    "predicted",
    "other",
]


def _default_accepted_interaction_types() -> list[InteractionType]:
    return ["physical", "direct", "genetic", "functional_association"]


class KnownInteractionsConfig(StrictModel):
    enabled: bool = False
    source: Literal["local_table"] = "local_table"
    local_table: Path | None = None
    minimum_supporting_records: int = Field(default=1, ge=1)
    minimum_direct_records: int = Field(default=1, ge=0)
    accepted_interaction_types: list[InteractionType] = Field(
        default_factory=_default_accepted_interaction_types
    )
    accepted_evidence_methods: list[str] = Field(default_factory=list)
    excluded_evidence_methods: list[str] = Field(default_factory=list)
    minimum_confidence: float | None = Field(default=None, ge=0.0, le=1.0)
    external_services: list[str] = Field(default_factory=list)


class ScoringWeightsConfig(StrictModel):
    genome_context: float = Field(default=1.0, ge=0.0)
    operon_proxy: float = Field(default=1.0, ge=0.0)
    domain_pair: float = Field(default=1.0, ge=0.0)
    functional_complementarity: float = Field(default=1.0, ge=0.0)
    localization: float = Field(default=0.5, ge=0.0)
    orthology: float = Field(default=0.75, ge=0.0)
    phylogenetic_profile: float = Field(default=1.0, ge=0.0)
    fusion: float = Field(default=1.5, ge=0.0)
    known_interactions: float = Field(default=1.5, ge=0.0)


class ScoringCategoryCapsConfig(StrictModel):
    genomic_context: float = Field(default=1.5, gt=0.0)
    functional_annotation: float = Field(default=1.5, gt=0.0)
    cellular_compatibility: float = Field(default=0.5, gt=0.0)
    evolutionary: float = Field(default=2.0, gt=0.0)
    direct_interaction: float = Field(default=2.0, gt=0.0)


class ScoringPenaltiesConfig(StrictModel):
    contradictory_evidence: float = Field(default=0.25, ge=0.0)
    ambiguous_mapping: float = Field(default=0.10, ge=0.0)


class ScoringConfig(StrictModel):
    enabled: bool = False
    rule_version: Literal["mvp1k-integrated-scoring-v1"] = "mvp1k-integrated-scoring-v1"
    output_scale: float = Field(default=100.0, gt=0.0)
    minimum_evidence_weight: float = Field(default=1.0, gt=0.0)
    minimum_evidence_categories: int = Field(default=2, ge=1)
    missing_policy: Literal["exclude_from_denominator"] = "exclude_from_denominator"
    tie_precision: int = Field(default=8, ge=0)
    weights: ScoringWeightsConfig = Field(default_factory=ScoringWeightsConfig)
    category_caps: ScoringCategoryCapsConfig = Field(default_factory=ScoringCategoryCapsConfig)
    penalties: ScoringPenaltiesConfig = Field(default_factory=ScoringPenaltiesConfig)


class StructurePredictionQueueConfig(StrictModel):
    enabled: bool = True
    maximum_entries: int = Field(default=20, ge=1)
    write_pair_fasta: bool = True
    automatic_structure_prediction: Literal[False] = False


class OutputConfig(StrictModel):
    directory: Path = Path("output")
    write_jsonl: bool = True
    write_tsv: bool = True
    write_excel: bool = True
    write_config_snapshot: bool = True


class CacheConfig(StrictModel):
    enabled: bool = True
    directory: Path = Path(".cache")


class LoggingConfig(StrictModel):
    level: Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"] = "INFO"
    directory: Path = Path("logs")


class PerformanceConfig(StrictModel):
    workers: int = Field(default=1, ge=1)


class AppConfig(StrictModel):
    project: ProjectConfig
    input: InputConfig
    query: QueryConfig
    candidate_generation: CandidateGenerationConfig
    gene_context: GeneContextConfig
    orthology: OrthologyConfig
    phylogenetic_profile: PhylogeneticProfileConfig = Field(
        default_factory=PhylogeneticProfileConfig
    )
    fusion: FusionConfig = Field(default_factory=FusionConfig)
    domains: DomainConfig
    functional_complementarity: FunctionalComplementarityConfig
    localization: LocalizationConfig
    known_interactions: KnownInteractionsConfig = Field(default_factory=KnownInteractionsConfig)
    scoring: ScoringConfig = Field(default_factory=ScoringConfig)
    evidence_tiers: EnabledConfig
    structure_prediction_queue: StructurePredictionQueueConfig
    output: OutputConfig
    cache: CacheConfig
    logging: LoggingConfig
    performance: PerformanceConfig


def _resolve(path: Path | None, base_directory: Path) -> Path | None:
    if path is None or path.is_absolute():
        return path
    return (base_directory / path).resolve()


def resolve_config_paths(config: AppConfig, base_directory: Path) -> AppConfig:
    """Return a copy whose paths are relative to the YAML file directory."""
    data = config.model_dump()
    input_data = data["input"]
    for key in ("proteome_fasta", "genome_gff", "annotation_table"):
        input_data[key] = _resolve(input_data[key], base_directory)
    input_data["reference_proteomes"] = [
        _resolve(path, base_directory) for path in input_data["reference_proteomes"]
    ]
    data["orthology"]["database"] = _resolve(
        data["orthology"]["database"],
        base_directory,
    )
    data["orthology"]["local_table"] = _resolve(
        data["orthology"]["local_table"],
        base_directory,
    )
    data["phylogenetic_profile"]["local_table"] = _resolve(
        data["phylogenetic_profile"]["local_table"],
        base_directory,
    )
    data["fusion"]["local_table"] = _resolve(
        data["fusion"]["local_table"],
        base_directory,
    )

    data["domains"]["local_table"] = _resolve(
        data["domains"]["local_table"],
        base_directory,
    )
    data["domains"]["rules_path"] = _resolve(
        data["domains"]["rules_path"],
        base_directory,
    )

    data["functional_complementarity"]["rules_path"] = _resolve(
        data["functional_complementarity"]["rules_path"],
        base_directory,
    )
    data["known_interactions"]["local_table"] = _resolve(
        data["known_interactions"]["local_table"], base_directory
    )
    for section in ("output", "cache", "logging"):
        data[section]["directory"] = _resolve(data[section]["directory"], base_directory)
    return AppConfig.model_validate(data)


def load_config(path: Path) -> AppConfig:
    """Load and strictly validate YAML without checking optional file existence."""
    config_path = path.expanduser().resolve()
    if not config_path.is_file():
        raise ConfigurationError(f"Configuration file not found: {config_path}")
    try:
        raw = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise ConfigurationError(f"Could not read YAML configuration: {exc}") from exc
    if not isinstance(raw, dict):
        raise ConfigurationError("Configuration root must be a YAML mapping.")
    try:
        config = AppConfig.model_validate(raw)
    except ValueError as exc:
        raise ConfigurationError(str(exc)) from exc
    return resolve_config_paths(config, config_path.parent)
