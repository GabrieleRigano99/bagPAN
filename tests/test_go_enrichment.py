from bagpan.go_enrichment import (
    compute_go_enrichment,
    orthogroup_go_associations,
    parse_obo,
    propagate_go_terms,
    run_go_enrichment,
)
from bagpan.orthogroups import classify_orthogroups

OBO_TEXT = """format-version: 1.2

[Term]
id: GO:0000003
name: root process

[Term]
id: GO:0000002
name: middle process
is_a: GO:0000003 ! root process

[Term]
id: GO:0000001
name: leaf process
is_a: GO:0000002 ! middle process
alt_id: GO:0000099

[Term]
id: GO:0000042
name: obsolete thing
is_obsolete: true

[Term]
id: GO:0000010
name: part-of example
relationship: part_of GO:0000003 ! root process
"""


def test_parse_obo_builds_parents_names_and_alts(tmp_path):
    obo_path = tmp_path / "go-basic.obo"
    obo_path.write_text(OBO_TEXT)
    dag = parse_obo(obo_path)

    assert dag.parents["GO:0000001"] == {"GO:0000002"}
    assert dag.parents["GO:0000002"] == {"GO:0000003"}
    assert dag.parents["GO:0000003"] == set()
    assert dag.parents["GO:0000010"] == {"GO:0000003"}  # part_of counted as a parent too
    assert dag.names["GO:0000001"] == "leaf process"
    assert dag.alt_to_primary["GO:0000099"] == "GO:0000001"
    assert "GO:0000042" not in dag.parents  # obsolete term dropped


def test_propagate_go_terms_expands_to_ancestors(tmp_path):
    obo_path = tmp_path / "go-basic.obo"
    obo_path.write_text(OBO_TEXT)
    dag = parse_obo(obo_path)

    propagated = propagate_go_terms({"protA": {"GO:0000001"}}, dag)
    assert propagated["protA"] == {"GO:0000001", "GO:0000002", "GO:0000003"}


def test_propagate_go_terms_resolves_alt_ids(tmp_path):
    obo_path = tmp_path / "go-basic.obo"
    obo_path.write_text(OBO_TEXT)
    dag = parse_obo(obo_path)

    propagated = propagate_go_terms({"protA": {"GO:0000099"}}, dag)
    assert propagated["protA"] == {"GO:0000001", "GO:0000002", "GO:0000003"}


def test_propagate_go_terms_none_dag_is_noop():
    terms = {"protA": {"GO:0000001"}}
    assert propagate_go_terms(terms, None) is terms


LOCUS_PREFIX_TO_SPECIES = {"W": "sw", "X": "sx"}
SPECIES = ["sw", "sx"]


def _five_core_five_singleton_scenario():
    clusters = {}
    for i in range(1, 6):
        clusters[f"CORE{i}"] = [f"W_{i}", f"X_{i}"]
    for i in range(1, 6):
        clusters[f"SINGLE{i}"] = [f"W_s{i}"]
    orthogroup_set = classify_orthogroups(clusters, LOCUS_PREFIX_TO_SPECIES, SPECIES)

    go_by_protein = {}
    for i in range(1, 6):
        go_by_protein[f"W_{i}"] = {"GO:0000001"}
        go_by_protein[f"X_{i}"] = {"GO:0000001"}
    return orthogroup_set, go_by_protein


def test_orthogroup_go_associations_respects_percent_threshold():
    orthogroup_set, go_by_protein = _five_core_five_singleton_scenario()
    associations = orthogroup_go_associations(orthogroup_set, go_by_protein, LOCUS_PREFIX_TO_SPECIES, percent=0.5)
    assert associations["CORE1"] == {"GO:0000001"}
    assert associations["SINGLE1"] == set()


def test_compute_go_enrichment_flags_perfectly_core_associated_term():
    orthogroup_set, go_by_protein = _five_core_five_singleton_scenario()
    associations = orthogroup_go_associations(orthogroup_set, go_by_protein, LOCUS_PREFIX_TO_SPECIES, percent=0.5)
    results = compute_go_enrichment(orthogroup_set, associations, min_population_count=1)

    core_row = next(r for r in results["Core"] if r[0] == "GO:0000001")
    go_id, study_hit, study_total, pop_hit, pop_total, p, q = core_row
    assert study_hit == 5 and study_total == 5
    assert pop_hit == 5 and pop_total == 10
    assert p < 0.05

    # a term absent everywhere in the singleton study set, present only elsewhere
    singleton_row = next(r for r in results["Singleton"] if r[0] == "GO:0000001")
    assert singleton_row[1] == 0  # study_hit


def test_compute_go_enrichment_min_population_count_filters_rare_terms():
    orthogroup_set, go_by_protein = _five_core_five_singleton_scenario()
    associations = orthogroup_go_associations(orthogroup_set, go_by_protein, LOCUS_PREFIX_TO_SPECIES, percent=0.5)
    results = compute_go_enrichment(orthogroup_set, associations, min_population_count=6)
    assert results["Core"] == []  # GO:0000001 only occurs in 5 orthogroups, below the cutoff


def test_run_go_enrichment_writes_one_file_per_class(tmp_path):
    orthogroup_set, go_by_protein = _five_core_five_singleton_scenario()
    run_go_enrichment(
        tmp_path, orthogroup_set, go_by_protein, LOCUS_PREFIX_TO_SPECIES, percent=0.5, alpha=0.05, min_population_count=1
    )
    go_dir = tmp_path / "go_enrichment"
    for cls in ("core", "accessory", "singleton"):
        assert (go_dir / f"go_enrichment_{cls}.tsv").exists()

    content = (go_dir / "go_enrichment_core.tsv").read_text()
    assert "GO:0000001" in content
    assert "Yes" in content  # enriched column should flag it given the strong signal