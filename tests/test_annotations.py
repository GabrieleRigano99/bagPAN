from pathlib import Path

from bagpan import annotations as ann
from bagpan.discover import discover_species

FIXTURES = Path(__file__).parent / "fixtures"


def test_parse_functional_annotation_tsv_basic_fields():
    sf = discover_species("species_a", FIXTURES / "species_a")
    genes = ann.parse_functional_annotation_tsv(sf.functional_annotation_tsv)

    assert set(genes) == {"SPA_000001", "SPA_000002", "SPA_000003"}
    gene = genes["SPA_000001"]
    assert gene.product == "Foo protein 1"
    assert gene.mrna_ids == {"SPA_000001-T1"}
    assert gene.pfam == {"PF00001"}
    assert gene.secreted is True
    assert gene.signalp is True
    assert gene.interpro == set()
    assert gene.go_terms == set()


def test_parse_functional_annotation_tsv_multivalue_fields():
    sf = discover_species("species_b", FIXTURES / "species_b")
    genes = ann.parse_functional_annotation_tsv(sf.functional_annotation_tsv)

    gene = genes["SPB_000001"]
    assert gene.interpro == {"IPR000001", "IPR001138"}
    assert gene.go_terms == {"GO:0003700"}
    assert gene.secreted is True
    assert gene.kegg_pathways == {"map00062", "map01100"}


def test_split_kegg_pathways_dedups_ko_and_map_prefixes():
    assert ann._split_kegg_pathways("ko00062|map00062") == {"map00062"}
    assert ann._split_kegg_pathways("ko00062|ko01100|map00062|map01100") == {"map00062", "map01100"}
    assert ann._split_kegg_pathways("") == set()
    # unrecognized format is kept as-is rather than silently dropped
    assert ann._split_kegg_pathways("weirdformat123") == {"weirdformat123"}


def test_parse_functional_annotation_tsv_transmembrane_and_effector():
    sf = discover_species("species_c", FIXTURES / "species_c")
    genes = ann.parse_functional_annotation_tsv(sf.functional_annotation_tsv)

    secreted_gene = genes["SPC_000001"]
    assert secreted_gene.secreted is False
    assert secreted_gene.pfam == {"PF00001"}

    tm_gene = genes["SPC_000002"]
    assert tm_gene.tm_phobius == 3
    assert tm_gene.tm_tmbed == 0
    assert tm_gene.effector_class == "Apoplastic effector"


def test_parse_functional_annotation_tsv_missing_file_returns_empty():
    assert ann.parse_functional_annotation_tsv(None) == {}


def test_parse_effectorp3_new_column_format(tmp_path):
    path = tmp_path / "effectorp3_output.txt"
    path.write_text(
        "# Identifier\tCytoplasmic effector\tApoplastic effector\tNon-effector\tPrediction\n"
        "PROT_1 gene=PROT_1 seq_id=ctg1 type=cds\tY (0.97)\t-\t-\tCytoplasmic effector\n"
        "PROT_2 gene=PROT_2 seq_id=ctg1 type=cds\t-\t-\tY (0.9)\tNon-effector\n"
    )
    effectors = ann.parse_effectorp3(path)
    assert set(effectors) == {"PROT_1"}
    assert abs(effectors["PROT_1"] - 0.97) < 1e-9


def test_parse_effectorp3_missing_file_returns_empty():
    assert ann.parse_effectorp3(None) == {}


def test_parse_annotation_stats(tmp_path):
    path = tmp_path / "annotation_stats.txt"
    path.write_text(
        "=== bagRNA Functional Annotation Summary ===\n\n"
        "Total genes                        17136  (100.0%)\n"
        "With product name                   4187  (24.4%)\n"
        "Secreted proteins                    512  (3.0%)\n"
    )
    stats = ann.parse_annotation_stats(path)
    assert stats["Total genes"] == (17136, "100.0%")
    assert stats["With product name"] == (4187, "24.4%")
    assert stats["Secreted proteins"] == (512, "3.0%")


def test_parse_annotation_stats_missing_file_returns_empty():
    assert ann.parse_annotation_stats(None) == {}
