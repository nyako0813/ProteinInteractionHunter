"""MVP-1B orchestration through candidate generation and observed gene context."""

from __future__ import annotations

import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path

import yaml
from pydantic import Field

from protein_interaction_hunter.adapters.local.annotation import LocalAnnotationTsvLoader
from protein_interaction_hunter.adapters.local.fasta import LocalFastaLoader
from protein_interaction_hunter.adapters.local.gff import LocalGff3Loader
from protein_interaction_hunter.application.candidates import generate_candidates
from protein_interaction_hunter.application.gene_context import (
    GENE_CONTEXT_RULE_VERSION,
    build_coordinate_index,
    calculate_gene_context,
)
from protein_interaction_hunter.application.identifiers import NORMALIZATION_RULE_VERSION
from protein_interaction_hunter.config import AppConfig, load_config
from protein_interaction_hunter.exceptions import InputValidationError
from protein_interaction_hunter.manifest import build_input_file_manifest, build_run_manifest
from protein_interaction_hunter.models.annotation import AnnotationRecord
from protein_interaction_hunter.models.base import StrictModel
from protein_interaction_hunter.models.enums import (
    CandidateDisposition,
    ContextCompleteness,
    EvidenceStatus,
    PredictedRelationshipType,
    RunStatus,
)
from protein_interaction_hunter.models.evidence import (
    CandidateEvidenceBundle,
    EvidenceProvenance,
    GenomeContextEvidence,
)
from protein_interaction_hunter.outputs.candidates import (
    CandidateTableTsvWriter,
    WarningSummaryTsvWriter,
)
from protein_interaction_hunter.outputs.excel import ExcelSchemaWriter
from protein_interaction_hunter.outputs.jsonl import (
    JsonlEvidenceBundleWriter,
    JsonRunManifestWriter,
)

_UNIMPLEMENTED_ENGINES = (
    "operon",
    "orthology",
    "phylogenetic_profile",
    "domains",
    "functional_complementarity",
    "localization",
    "fusion",
    "known_interactions",
    "scoring",
    "evidence_tiers",
)


class CandidateGenerationSummary(StrictModel):
    query_count: int = Field(ge=0)
    protein_count: int = Field(ge=0)
    pair_count: int = Field(ge=0)
    included_count: int = Field(ge=0)
    flagged_count: int = Field(ge=0)
    excluded_count: int = Field(ge=0)
    duplicate_group_count: int = Field(ge=0)
    fragment_candidate_count: int = Field(ge=0)
    hypothetical_protein_count: int = Field(ge=0)
    missing_coordinate_count: int = Field(ge=0)
    missing_annotation_count: int = Field(ge=0)
    ambiguous_mapping_count: int = Field(ge=0)
    same_contig_pair_count: int = Field(default=0, ge=0)
    different_contig_pair_count: int = Field(default=0, ge=0)
    overlapping_pair_count: int = Field(default=0, ge=0)
    missing_context_pair_count: int = Field(default=0, ge=0)
    ambiguous_context_pair_count: int = Field(default=0, ge=0)
    neighborhood_pair_count: int = Field(default=0, ge=0)
    incomplete_context_pair_count: int = Field(default=0, ge=0)
    output_path: Path


class PipelineResult(StrictModel):
    summary: CandidateGenerationSummary
    bundles: list[CandidateEvidenceBundle]
    evidence_path: Path
    candidate_table_path: Path
    manifest_path: Path
    config_snapshot_path: Path
    warning_summary_path: Path
    excel_path: Path | None = None


def _git_commit(repository: Path) -> str | None:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repository,
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        return None
    return result.stdout.strip() or None


def _policy_settings(config: AppConfig) -> dict[str, str | int | bool]:
    policy = config.candidate_generation
    return {
        "include_hypothetical_proteins": policy.include_hypothetical_proteins,
        "include_missing_coordinates": policy.include_missing_coordinates,
        "minimum_length_aa": policy.minimum_length_aa,
        "duplicate_sequence_policy": policy.duplicate_sequence_policy,
        "fragment_policy": policy.fragment_policy,
        "self_candidate_policy": policy.self_candidate_policy,
        "gene_context_enabled": config.gene_context.enabled,
        "neighborhood_gene_count": config.gene_context.neighborhood_gene_count,
        "require_query_coordinates": config.gene_context.require_query_coordinates,
    }


def _write_snapshot(config: AppConfig, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        yaml.safe_dump(config.model_dump(mode="json"), sort_keys=True, allow_unicode=True),
        encoding="utf-8",
        newline="\n",
    )
    return path


def _excel_context_rows(
    run_id: str,
    contexts: dict[tuple[str, str], GenomeContextEvidence],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for (query_id, candidate_id), context in sorted(contexts.items()):
        rows.append(
            {
                "Run_ID": run_id,
                "Query_ID": query_id,
                "Candidate_ID": candidate_id,
                "Same_Contig": context.same_contig,
                "Query_Contig": context.query_contig,
                "Candidate_Contig": context.candidate_contig,
                "Query_Start": context.query_start,
                "Query_End": context.query_end,
                "Query_Strand": context.query_strand,
                "Candidate_Start": context.candidate_start,
                "Candidate_End": context.candidate_end,
                "Candidate_Strand": context.candidate_strand,
                "Strand_Relationship": context.strand_relationship,
                "Relative_Position": context.relative_position,
                "Coordinate_Position": context.coordinate_position,
                "Distance_BP": context.distance_bp,
                "Overlap_BP": context.overlap_bp,
                "Intervening_Gene_Count": context.intervening_gene_count,
                "Intervening_Feature_Count": context.intervening_feature_count,
                "Feature_Index_Delta": context.feature_index_delta,
                "Within_Neighborhood_Window": context.within_neighborhood_window,
                "Context_Completeness": context.context_completeness,
                "Status": context.status,
                "Warnings": "|".join(context.warnings),
                "Source": GENE_CONTEXT_RULE_VERSION,
                "Provenance": "|".join(
                    f"{item.source_name}:{item.source_version or ''}:{item.method or ''}"
                    for item in context.provenance
                ),
            }
        )
    return rows


class InteractionCandidatePipeline:
    """Generate auditable candidates and coordinate-derived context without scoring."""

    def run(
        self, config_path: Path | None = None, command_line: list[str] | None = None
    ) -> PipelineResult:
        if config_path is None:
            raise NotImplementedError(
                "Candidate ranking is not implemented in MVP-0; a config path is required."
            )
        resolved_config_path = config_path.expanduser().resolve()
        config = load_config(resolved_config_path)
        proteins = LocalFastaLoader().load(config.input.proteome_fasta)
        document = LocalGff3Loader().load_document(config.input.genome_gff)
        coordinates = document.features
        annotations: list[AnnotationRecord] = []
        if config.input.annotation_table is not None:
            annotations = LocalAnnotationTsvLoader().load(config.input.annotation_table)
        generated = generate_candidates(
            proteins=proteins,
            coordinates=coordinates,
            annotations=annotations,
            query_ids=config.query.protein_ids,
            policy=config.candidate_generation,
        )
        canonical_query_ids = {query.query_id: query.protein_id for query in generated.queries}
        contexts: dict[tuple[str, str], GenomeContextEvidence] = {}
        if config.gene_context.enabled:
            coordinate_index = build_coordinate_index(proteins, document)
            if config.gene_context.require_query_coordinates:
                for query in generated.queries:
                    if query.protein_id in coordinate_index.ambiguous_proteins:
                        raise InputValidationError(
                            f"Query coordinate is ambiguous: {query.protein_id}"
                        )
                    if query.protein_id not in coordinate_index.by_protein:
                        raise InputValidationError(
                            f"Query coordinate is missing: {query.protein_id}"
                        )
            for candidate in generated.candidates:
                contexts[(candidate.query_id, candidate.protein_id)] = calculate_gene_context(
                    canonical_query_ids[candidate.query_id],
                    candidate.protein_id,
                    coordinate_index,
                    config.gene_context.neighborhood_gene_count,
                )
        config_hash = build_input_file_manifest(
            "config", resolved_config_path, required=True
        ).sha256
        assert config_hash is not None
        run_id = f"{config.project.run_name}-{config_hash[:12]}"
        output_path = config.output.directory
        output_path.mkdir(parents=True, exist_ok=True)
        snapshot_path = _write_snapshot(config, output_path / "config.snapshot.yaml")
        policies = _policy_settings(config)
        provenance = EvidenceProvenance(
            source_name="candidate_generation",
            source_version=NORMALIZATION_RULE_VERSION,
            method="deterministic_all_proteome_enumeration",
            metadata={"policy_settings": policies},
        )
        bundles: list[CandidateEvidenceBundle] = []
        for candidate in generated.candidates:
            context = contexts.get((candidate.query_id, candidate.protein_id))
            statuses = {engine: EvidenceStatus.NOT_RUN for engine in _UNIMPLEMENTED_ENGINES}
            statuses["gene_context"] = context.status if context else EvidenceStatus.NOT_RUN
            context_warnings = context.warnings if context else []
            bundles.append(
                CandidateEvidenceBundle(
                    run_id=run_id,
                    query_id=candidate.query_id,
                    candidate_id=candidate.protein_id,
                    candidate=candidate,
                    candidate_disposition=candidate.disposition,
                    predicted_relationship_type=PredictedRelationshipType.INSUFFICIENT_EVIDENCE,
                    genome_context=[context] if context else [],
                    engine_statuses=statuses,
                    provenance=[provenance],
                    policy_settings=policies,
                    warnings=sorted(set(candidate.warnings + context_warnings)),
                )
            )
        evidence_path = JsonlEvidenceBundleWriter().write(
            bundles, output_path / "candidate_evidence_bundle.jsonl"
        )
        candidate_table_path = CandidateTableTsvWriter().write(
            run_id,
            generated.candidates,
            output_path / "candidate_table.tsv",
            contexts=contexts,
        )
        all_warnings = [warning for bundle in bundles for warning in bundle.warnings]
        warning_summary_path = WarningSummaryTsvWriter().write(
            all_warnings, output_path / "warning_summary.tsv"
        )
        excel_path = None
        if config.output.write_excel:
            excel_path = ExcelSchemaWriter().write(
                output_path / "ProteinInteractionHunter.xlsx",
                rows_by_sheet={"Gene_Context": _excel_context_rows(run_id, contexts)},
            )
        input_files = [
            build_input_file_manifest("proteome_fasta", config.input.proteome_fasta, required=True),
            build_input_file_manifest("genome_gff", config.input.genome_gff, required=True),
        ]
        if config.input.annotation_table is not None:
            input_files.append(
                build_input_file_manifest(
                    "annotation_table", config.input.annotation_table, required=False
                )
            )
        manifest = build_run_manifest(
            run_id=run_id,
            run_name=config.project.run_name,
            config_path=resolved_config_path,
            config_snapshot_path=snapshot_path,
            input_files=input_files,
            random_seed=config.project.random_seed,
            command_line=command_line or sys.argv,
        )
        manifest.status = RunStatus.COMPLETED
        manifest.completed_at = datetime.now(UTC)
        manifest.git_commit = _git_commit(Path(__file__).resolve().parents[3])
        manifest.normalization_rule_version = NORMALIZATION_RULE_VERSION
        manifest.gene_context_rule_version = (
            GENE_CONTEXT_RULE_VERSION if config.gene_context.enabled else None
        )
        manifest.policy_settings = policies
        manifest.warnings = sorted(set(all_warnings))
        manifest.parser_warnings = sorted(
            set(document.warnings)
            | {warning for protein in proteins for warning in protein.warnings}
            | {warning for coordinate in coordinates for warning in coordinate.warnings}
            | {warning for annotation in annotations for warning in annotation.warnings}
        )
        manifest.incomplete_evidence_flags = [
            f"{engine}_not_run" for engine in _UNIMPLEMENTED_ENGINES
        ]
        if not config.gene_context.enabled:
            manifest.incomplete_evidence_flags.append("gene_context_not_run")
        manifest_path = JsonRunManifestWriter().write(manifest, output_path / "run_manifest.json")
        disposition_counts = {
            disposition: sum(
                candidate.disposition is disposition for candidate in generated.candidates
            )
            for disposition in CandidateDisposition
        }
        context_values = list(contexts.values())
        summary = CandidateGenerationSummary(
            query_count=len(generated.queries),
            protein_count=len(proteins),
            pair_count=len(generated.candidates),
            included_count=disposition_counts[CandidateDisposition.INCLUDED],
            flagged_count=disposition_counts[CandidateDisposition.FLAGGED],
            excluded_count=disposition_counts[CandidateDisposition.EXCLUDED],
            duplicate_group_count=len(generated.duplicate_groups),
            fragment_candidate_count=sum(
                candidate.is_fragment_candidate for candidate in generated.candidates
            ),
            hypothetical_protein_count=len(
                {
                    candidate.protein_id
                    for candidate in generated.candidates
                    if candidate.is_hypothetical
                }
            ),
            missing_coordinate_count=len(
                {
                    candidate.protein_id
                    for candidate in generated.candidates
                    if not candidate.has_coordinate
                }
            ),
            missing_annotation_count=len(
                {
                    candidate.protein_id
                    for candidate in generated.candidates
                    if not candidate.has_annotation
                }
            ),
            ambiguous_mapping_count=generated.ambiguous_mapping_count,
            same_contig_pair_count=sum(item.same_contig is True for item in context_values),
            different_contig_pair_count=sum(
                item.status is EvidenceStatus.NOT_APPLICABLE for item in context_values
            ),
            overlapping_pair_count=sum(
                item.status is EvidenceStatus.AVAILABLE and bool(item.overlap_bp)
                for item in context_values
            ),
            missing_context_pair_count=sum(
                item.status is EvidenceStatus.MISSING for item in context_values
            ),
            ambiguous_context_pair_count=sum(
                item.status is EvidenceStatus.FAILED for item in context_values
            ),
            neighborhood_pair_count=sum(
                item.within_neighborhood_window is True for item in context_values
            ),
            incomplete_context_pair_count=sum(
                item.status is EvidenceStatus.AVAILABLE
                and item.context_completeness is not ContextCompleteness.COMPLETE
                for item in context_values
            ),
            output_path=output_path.resolve(),
        )
        return PipelineResult(
            summary=summary,
            bundles=bundles,
            evidence_path=evidence_path,
            candidate_table_path=candidate_table_path,
            manifest_path=manifest_path,
            config_snapshot_path=snapshot_path,
            warning_summary_path=warning_summary_path,
            excel_path=excel_path,
        )
