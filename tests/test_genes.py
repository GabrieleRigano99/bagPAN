from pathlib import Path

from bagpan.discover import discover_species
from bagpan.genes import build_gene_annotation, _from_gff3, _from_protein_ids_fallback

FIXTURES = Path(__file__).parent / "fixtures"


def test_build_gene_annotation_from_real_gff3():
    sf = discover_species("species_a", FIXTURES / "species_a")
    assert sf.gff3 is not None
    ga = build_gene_annotation(sf)
    assert ga.from_gff3 is True
    assert ga.gene_of["SPA_000001-T1"] == "SPA_000001"
    assert ga.representative_of_gene["SPA_000001"] == "SPA_000001-T1"
    assert ga.gene_position["SPA_000001"] == ("ctg1", 100, "+")
    assert ga.gene_position["SPA_000002"] == ("ctg1", 500, "+")
    assert ga.transcript_contig["SPA_000001-T1"] == "ctg1"
    assert ga.transcript_strand["SPA_000001-T1"] == "+"
    assert ga.transcript_cds_intervals["SPA_000001-T1"] == [(100, 400)]


def test_build_gene_annotation_captures_multi_exon_cds_intervals(tmp_path):
    from bagpan.genes import _from_gff3

    gff3 = tmp_path / "toy.gff3"
    gff3.write_text(
        "##gff-version 3\n"
        "ctg1\tx\tgene\t1\t1000\t.\t-\t.\tID=GENE1;\n"
        "ctg1\tx\tmRNA\t1\t1000\t.\t-\t.\tID=GENE1-T1;Parent=GENE1;\n"
        "ctg1\tx\tCDS\t600\t1000\t.\t-\t0\tID=GENE1-T1.cds;Parent=GENE1-T1;\n"
        "ctg1\tx\tCDS\t1\t400\t.\t-\t0\tID=GENE1-T1.cds;Parent=GENE1-T1;\n"
    )
    ga = _from_gff3(gff3)
    assert ga.transcript_strand["GENE1-T1"] == "-"
    assert ga.transcript_contig["GENE1-T1"] == "ctg1"
    assert sorted(ga.transcript_cds_intervals["GENE1-T1"]) == [(1, 400), (600, 1000)]


def test_representative_transcript_picks_longest_cds(tmp_path):
    gff3 = tmp_path / "toy.gff3"
    gff3.write_text(
        "##gff-version 3\n"
        "ctg1\tx\tgene\t1\t1000\t.\t+\t.\tID=GENE1;\n"
        "ctg1\tx\tmRNA\t1\t500\t.\t+\t.\tID=GENE1-T1;Parent=GENE1;\n"
        "ctg1\tx\tCDS\t1\t500\t.\t+\t0\tID=GENE1-T1.cds;Parent=GENE1-T1;\n"
        "ctg1\tx\tmRNA\t1\t1000\t.\t+\t.\tID=GENE1-T2;Parent=GENE1;\n"
        "ctg1\tx\tCDS\t1\t400\t.\t+\t0\tID=GENE1-T2.cds;Parent=GENE1-T2;\n"
        "ctg1\tx\tCDS\t600\t1000\t.\t+\t0\tID=GENE1-T2.cds;Parent=GENE1-T2;\n"
    )
    ga = _from_gff3(gff3)
    # T1 CDS length = 500; T2 CDS length = 400 + 401 = 801 -> T2 should win
    assert ga.representative_of_gene["GENE1"] == "GENE1-T2"


def test_representative_transcript_ties_broken_deterministically(tmp_path):
    gff3 = tmp_path / "toy.gff3"
    gff3.write_text(
        "##gff-version 3\n"
        "ctg1\tx\tgene\t1\t1000\t.\t+\t.\tID=GENE1;\n"
        "ctg1\tx\tmRNA\t1\t500\t.\t+\t.\tID=GENE1-T2;Parent=GENE1;\n"
        "ctg1\tx\tCDS\t1\t500\t.\t+\t0\tID=GENE1-T2.cds;Parent=GENE1-T2;\n"
        "ctg1\tx\tmRNA\t1\t500\t.\t+\t.\tID=GENE1-T1;Parent=GENE1;\n"
        "ctg1\tx\tCDS\t1\t500\t.\t+\t0\tID=GENE1-T1.cds;Parent=GENE1-T1;\n"
    )
    ga = _from_gff3(gff3)
    assert ga.representative_of_gene["GENE1"] == "GENE1-T1"


def test_fallback_strips_transcript_suffix_when_no_gff3():
    ga = _from_protein_ids_fallback(["FOO_0001-T1", "FOO_0001-T2", "FOO_0002-T1"])
    assert ga.from_gff3 is False
    assert ga.gene_of["FOO_0001-T1"] == "FOO_0001"
    assert ga.gene_of["FOO_0001-T2"] == "FOO_0001"
    assert ga.representative_of_gene["FOO_0001"] == "FOO_0001-T1"
    assert ga.representative_of_gene["FOO_0002"] == "FOO_0002-T1"


def test_build_gene_annotation_falls_back_without_gff3(tmp_path, capsys):
    from bagpan.discover import SpeciesFiles

    proteome = tmp_path / "toy.proteins.fa"
    proteome.write_text(">FOO_0001-T1 desc\nMSEQ\n>FOO_0001-T2 desc\nMSEQ\n")
    sf = SpeciesFiles(
        name="toy", root=tmp_path, proteome=proteome, locus_prefix="FOO",
        functional_annotation_tsv=tmp_path / "functional_annotation.tsv", gff3=None,
    )
    ga = build_gene_annotation(sf)
    assert ga.from_gff3 is False
    assert ga.gene_of["FOO_0001-T1"] == "FOO_0001"
    assert "WARNING" in capsys.readouterr().out