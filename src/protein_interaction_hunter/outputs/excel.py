"""Minimal schema workbook writer; it performs no candidate analysis."""

from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

from openpyxl import Workbook

EXCEL_SHEETS: dict[str, tuple[str, ...]] = {
    "Run_Summary": ("Run_ID", "Run_Name", "Status", "Warning_Count"),
    "Input_Manifest": (
        "Logical_Name",
        "Path",
        "SHA256",
        "Size_Bytes",
        "Modified_Time",
        "Exists",
        "Required",
    ),
    "Query_Proteins": ("Query_ID", "Protein_ID", "Resolution_Method", "Warnings"),
    "Candidate_Ranking": (
        "Rank",
        "Query_ID",
        "Candidate_ID",
        "Disposition",
        "Relationship_Type",
        "Evidence_Tier",
        "Total_Ranking_Score",
    ),
    "Evidence_Detail": (
        "Query_ID",
        "Candidate_ID",
        "Category",
        "Status",
        "Origin",
        "Quality",
        "Source",
    ),
    "Gene_Context": (
        "Run_ID",
        "Query_ID",
        "Candidate_ID",
        "Same_Contig",
        "Query_Contig",
        "Candidate_Contig",
        "Query_Start",
        "Query_End",
        "Query_Strand",
        "Candidate_Start",
        "Candidate_End",
        "Candidate_Strand",
        "Strand_Relationship",
        "Relative_Position",
        "Coordinate_Position",
        "Distance_BP",
        "Overlap_BP",
        "Intervening_Gene_Count",
        "Intervening_Feature_Count",
        "Feature_Index_Delta",
        "Within_Neighborhood_Window",
        "Context_Completeness",
        "Status",
        "Warnings",
        "Source",
        "Provenance",
    ),
    "Operon_Proxy": (
        "Run_ID",
        "Query_ID",
        "Candidate_ID",
        "Status",
        "Proxy_Status",
        "Same_Contig",
        "Same_Strand",
        "Is_Adjacent",
        "Intergenic_Distance_BP",
        "Overlap_BP",
        "Intervening_Gene_Count",
        "Transcriptional_Order",
        "Maximum_Intergenic_Distance_BP",
        "Passes_Distance_Threshold",
        "Supporting_Conditions",
        "Conflicting_Conditions",
        "Rule_Version",
        "Rule_ID",
        "Warnings",
        "Provenance",
    ),
    "Contradictions": (
        "Query_ID",
        "Candidate_ID",
        "Type",
        "Severity",
        "Penalty",
        "Explanation",
    ),
    "Structure_Prediction_Queue": (
        "Rank",
        "Query_ID",
        "Candidate_ID",
        "Manual_Structure_Priority",
        "Manual_Review_Status",
    ),
    "Warnings": ("Stage", "Code", "Severity", "Entity_ID", "Message"),
    "Provenance": ("Component", "Version", "Source", "Retrieved_At", "Checksum"),
}


class ExcelSchemaWriter:
    def write(
        self,
        path: Path,
        rows_by_sheet: Mapping[str, Iterable[Mapping[str, Any]]] | None = None,
    ) -> Path:
        output_path = path.expanduser().resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        workbook = Workbook()
        active_sheet = workbook.active
        if active_sheet is not None:
            workbook.remove(active_sheet)
        supplied = rows_by_sheet or {}
        for sheet_name, headers in EXCEL_SHEETS.items():
            worksheet = workbook.create_sheet(sheet_name)
            worksheet.append(list(headers))
            worksheet.freeze_panes = "A2"
            for row in supplied.get(sheet_name, []):
                worksheet.append([row.get(header) for header in headers])
        workbook.save(output_path)
        return output_path
