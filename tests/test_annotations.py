from pathlib import Path

from bagpan import annotations as ann
from bagpan.discover import discover_species

FIXTURES = Path(__file__).parent / "fixtures"


def test_parse_generic_annotations_pfam_and_products():
    sf = discover_species("species_a", FIXTURES / "species_a")
    pfam = ann.parse_generic_annotations(sf.pfam)
    assert pfam["SPA_000001-T1"] == ["PF00001"]

    products = ann.parse_generic_annotations(sf.genes_products)
    assert products["SPA_000001-T1"] == ["Foo protein 1"]


def test_parse_go_skips_non_go_lines_and_extracts_term():
    sf = discover_species("species_b", FIXTURES / "species_b")
    go = ann.parse_go(sf.iprscan)
    assert go["SPB_000001-T1"] == {"GO:0003700"}
    # db_xref line for the same protein must not leak into the GO set
    generic = ann.parse_generic_annotations(sf.iprscan)
    assert generic["SPB_000001-T1"] == ["IPR000001"]


def test_parse_phobius_secreted_vs_transmembrane():
    sf = discover_species("species_c", FIXTURES / "species_c")
    secreted, tm = ann.parse_phobius(sf.phobius)
    assert secreted == set()
    assert tm == {"SPC_000002-T1": 3}


def test_parse_signalp_secreted_ids():
    sf = discover_species("species_a", FIXTURES / "species_a")
    secreted = ann.parse_signalp(sf.signalp)
    assert secreted == {"SPA_000001-T1", "SPA_000002-T1"}


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


def test_parse_antismash_cleans_duplicated_prefix_and_cds_suffix():
    sf = discover_species("species_a", FIXTURES / "species_a")
    cluster_hits, smcog_hits = ann.parse_antismash([sf.antismash, sf.antismash_clusters])
    # the fixture's SMCOG line is 'SPA_SPA_000001-T1.cds' - must normalize to
    # the same id as the cluster-hit line 'SPA_000001-T1'
    assert smcog_hits == {"SPA_000001-T1"}
    assert cluster_hits == {"SPA_000001-T1"}


def test_clean_antismash_id_leaves_normal_ids_untouched():
    assert ann._clean_antismash_id("SPA_000002-T1") == "SPA_000002-T1"
    assert ann._clean_antismash_id("SPA_SPA_000005-T1.cds") == "SPA_000005-T1"