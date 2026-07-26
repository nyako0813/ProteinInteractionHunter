import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from scripts.build_ncbi_annotation_table import (  # noqa: E402
    build_annotation_tables,
    clean_fasta_description,
    write_build_outputs,
)


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="\n")
    return path


def test_one_to_one_query_mapping_and_missing_product(tmp_path: Path) -> None:
    fasta = _write(
        tmp_path / "proteins.faa",
        ">Q alpha hydrolase [Test archaeon]\nMSTKAA\n>M fallback product [Test archaeon]\nMSTKCC\n",
    )
    gff = _write(
        tmp_path / "genome.gff",
        "##gff-version 3\n"
        "##sequence-region c 1 1000\n"
        "c\tRefSeq\tgene\t10\t30\t.\t+\t.\t"
        "ID=gene-Q;gene=hydA;locus_tag=LTQ;old_locus_tag=OLDQ1%2COLDQ2;"
        "Dbxref=GeneID:1\n"
        "c\tRefSeq\tCDS\t10\t30\t.\t+\t0\t"
        "ID=cds-Q;Parent=gene-Q;protein_id=Q;locus_tag=LTQ;"
        "product=alpha%20hydrolase\n"
        "c\tRefSeq\tgene\t100\t120\t.\t+\t.\tID=gene-M;locus_tag=LTM\n"
        "c\tRefSeq\tCDS\t100\t120\t.\t+\t0\t"
        "ID=cds-M;Parent=gene-M;protein_id=M;locus_tag=LTM\n",
    )

    result = build_annotation_tables(
        fasta_path=fasta,
        gff_path=gff,
        annotation_source="NCBI Test",
        query_id="Q",
    )

    assert [row.protein_id for row in result.annotations] == ["M", "Q"]
    query = result.annotations[1]
    assert (query.gene_name, query.locus_tag, query.product) == (
        "hydA",
        "LTQ",
        "alpha hydrolase",
    )
    assert query.functional_category == query.localization_annotation == ""
    assert query.transmembrane_annotation == query.annotation_confidence == ""
    missing = result.annotations[0]
    assert missing.product == "fallback product"
    query_audit = next(row for row in result.audit_rows if row.protein_id == "Q")
    assert query_audit.old_locus_tags == "OLDQ1|OLDQ2"
    assert query_audit.ncbi_gene_id == "1"
    assert query_audit.product_match == "true"
    assert result.summary.query_covered is True
    assert result.summary.missing_product == 0


def test_fasta_only_and_gff_only_are_audited_without_unknown_formal_row(
    tmp_path: Path,
) -> None:
    fasta = _write(tmp_path / "proteins.faa", ">F fasta only [Taxon]\nMSTKAA\n")
    gff = _write(
        tmp_path / "genome.gff",
        "##gff-version 3\n"
        "c\tRefSeq\tCDS\t10\t30\t.\t+\t0\t"
        "ID=cds-G;protein_id=G;product=gff%20only\n",
    )

    result = build_annotation_tables(fasta_path=fasta, gff_path=gff, annotation_source="NCBI Test")

    assert [row.protein_id for row in result.annotations] == ["F"]
    assert result.annotations[0].product == "fasta only"
    assert {row.mapping_status for row in result.audit_rows} == {
        "fasta_only",
        "gff_only",
    }
    assert result.summary.fasta_only == result.summary.gff_only == 1


def test_multi_locus_accession_is_excluded_and_split_cds_is_merged(
    tmp_path: Path,
) -> None:
    fasta = _write(
        tmp_path / "proteins.faa",
        ">A repeated [Taxon]\nMSTKAA\n>S split product [Taxon]\nMSTKCC\n",
    )
    gff = _write(
        tmp_path / "genome.gff",
        "##gff-version 3\n"
        "c\tRefSeq\tCDS\t10\t30\t.\t+\t0\t"
        "ID=a1;Parent=gene-a1;protein_id=A;product=repeated\n"
        "c\tRefSeq\tCDS\t100\t120\t.\t+\t0\t"
        "ID=a2;Parent=gene-a2;protein_id=A;product=repeated\n"
        "c\tRefSeq\tCDS\t200\t210\t.\t-\t0\t"
        "ID=s1;Parent=gene-s;protein_id=S;product=split%20product\n"
        "c\tRefSeq\tCDS\t220\t230\t.\t-\t2\t"
        "ID=s2;Parent=gene-s;protein_id=S;product=split%20product\n",
    )

    result = build_annotation_tables(fasta_path=fasta, gff_path=gff, annotation_source="NCBI Test")

    assert [row.protein_id for row in result.annotations] == ["S"]
    ambiguous = [row for row in result.audit_rows if row.protein_id == "A"]
    assert len(ambiguous) == 2
    assert all(row.mapping_status == "ambiguous_multi_locus" for row in ambiguous)
    split = next(row for row in result.audit_rows if row.protein_id == "S")
    assert (split.start, split.end, split.cds_record_count) == (200, 230, 2)
    assert result.summary.multi_locus_ambiguous == 1


def test_duplicate_cds_is_deduplicated_semantically_and_reported(
    tmp_path: Path,
) -> None:
    fasta = _write(tmp_path / "proteins.faa", ">D duplicate [Taxon]\nMSTKAA\n")
    line = "c\tRefSeq\tCDS\t10\t30\t.\t+\t0\tID=d;Parent=gene-d;protein_id=D;product=duplicate\n"
    gff = _write(tmp_path / "genome.gff", "##gff-version 3\n" + line + line)

    result = build_annotation_tables(fasta_path=fasta, gff_path=gff, annotation_source="NCBI Test")

    assert [row.protein_id for row in result.annotations] == ["D"]
    assert result.audit_rows[0].duplicate_cds_count == 1
    assert "duplicate_cds_record" in result.audit_rows[0].mapping_warning
    assert result.summary.duplicate_cds_records == 1
    assert result.summary.duplicate_annotation_rows == 0


def test_output_is_deterministic_lf_and_input_order_independent(
    tmp_path: Path,
) -> None:
    fasta = _write(
        tmp_path / "proteins.faa",
        ">B beta [Taxon]\nMSTKAA\n>A alpha [Taxon]\nMSTKCC\n",
    )
    header = "##gff-version 3\n"
    lines = [
        "c\tRefSeq\tCDS\t100\t120\t.\t+\t0\tID=b;protein_id=B;locus_tag=LTB;product=beta\n",
        "c\tRefSeq\tCDS\t10\t30\t.\t+\t0\tID=a;protein_id=A;locus_tag=LTA;product=alpha\n",
    ]
    first_gff = _write(tmp_path / "first.gff", header + "".join(lines))
    second_gff = _write(tmp_path / "second.gff", header + "".join(reversed(lines)))
    first = build_annotation_tables(
        fasta_path=fasta, gff_path=first_gff, annotation_source="NCBI Test"
    )
    second = build_annotation_tables(
        fasta_path=fasta, gff_path=second_gff, annotation_source="NCBI Test"
    )
    assert first == second
    assert clean_fasta_description("alpha [Taxon]") == "alpha"
    assert clean_fasta_description("alpha [gene=x]") == "alpha [gene=x]"

    annotation = tmp_path / "annotation.tsv"
    audit = tmp_path / "audit.tsv"
    coverage = tmp_path / "coverage.tsv"
    write_build_outputs(
        first,
        annotation_output=annotation,
        audit_output=audit,
        coverage_output=coverage,
    )
    assert b"\r" not in annotation.read_bytes()
    assert b"\r" not in audit.read_bytes()
    assert b"\r" not in coverage.read_bytes()
    assert annotation.read_text(encoding="utf-8").splitlines()[1].startswith("A\t")
