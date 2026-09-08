from bagpan.dnds_runner import extract_all_cds, run_pairwise_dnds, select_single_copy_pairs
from bagpan.genes import GeneAnnotation
from bagpan.orthogroups import classify_orthogroups

LOCUS_PREFIX_TO_SPECIES = {"X": "sx", "Y": "sy", "Z": "sz"}
SPECIES = ["sx", "sy", "sz"]

CLUSTERS = {
    "OG_SC": ["X_1", "Y_1", "Z_1"],  # single-copy in all 3 - usable
    "OG_PARALOG": ["X_2", "X_3", "Y_2"],  # 2 copies in sx - not usable
    "OG_PAIR_ONLY": ["X_4", "Y_4"],  # single-copy, only 2 species present
    "OG_SINGLETON": ["X_5"],  # only 1 species - not usable (needs >=2)
}


def test_select_single_copy_pairs():
    orthogroup_set = classify_orthogroups(CLUSTERS, LOCUS_PREFIX_TO_SPECIES, SPECIES)
    result = select_single_copy_pairs(orthogroup_set, LOCUS_PREFIX_TO_SPECIES)
    assert set(result) == {"OG_SC", "OG_PAIR_ONLY"}
    assert result["OG_SC"] == {"sx": "X_1", "sy": "Y_1", "sz": "Z_1"}
    assert result["OG_PAIR_ONLY"] == {"sx": "X_4", "sy": "Y_4"}


def test_extract_all_cds_skips_missing_contigs_and_no_gff3_species():
    ga = GeneAnnotation(
        gene_of={"g1-T1": "g1", "g2-T1": "g2"},
        representative_of_gene={"g1": "g1-T1", "g2": "g2-T1"},
        gene_position={},
        from_gff3=True,
        transcript_contig={"g1-T1": "ctg1", "g2-T1": "ctg_missing"},
        transcript_strand={"g1-T1": "+", "g2-T1": "+"},
        transcript_cds_intervals={"g1-T1": [(1, 6)], "g2-T1": [(1, 6)]},
    )
    genome = {"ctg1": "ATGCCCTAA"}
    result = extract_all_cds(ga, genome)
    assert result == {"g1-T1": "ATGCCC"}  # g2-T1's contig isn't in the genome dict


def test_extract_all_cds_empty_for_fallback_gene_annotation():
    fallback_ga = GeneAnnotation(
        gene_of={"g1-T1": "g1"}, representative_of_gene={"g1": "g1-T1"},
        gene_position={}, from_gff3=False,
    )
    assert extract_all_cds(fallback_ga, {"ctg1": "ACGT"}) == {}


def test_run_pairwise_dnds_end_to_end_on_identical_sequences():
    orthogroup_set = classify_orthogroups(CLUSTERS, LOCUS_PREFIX_TO_SPECIES, SPECIES)
    cds = "ATG" + "GCT" * 20 + "TAA"
    cds_by_protein = {"X_1": cds, "Y_1": cds, "Z_1": cds, "X_4": cds, "Y_4": cds}

    rows = run_pairwise_dnds(orthogroup_set, LOCUS_PREFIX_TO_SPECIES, cds_by_protein, max_orthogroups=200)

    og_sc_rows = [r for r in rows if r["orthogroup"] == "OG_SC"]
    assert len(og_sc_rows) == 3  # C(3,2) species pairs
    for row in og_sc_rows:
        assert row["dS"] == 0.0 and row["dN"] == 0.0

    og_pair_rows = [r for r in rows if r["orthogroup"] == "OG_PAIR_ONLY"]
    assert len(og_pair_rows) == 1

    # paralog/singleton orthogroups never appear
    assert not any(r["orthogroup"] in ("OG_PARALOG", "OG_SINGLETON") for r in rows)


def test_run_pairwise_dnds_skips_pairs_with_missing_cds():
    orthogroup_set = classify_orthogroups(CLUSTERS, LOCUS_PREFIX_TO_SPECIES, SPECIES)
    cds = "ATG" + "GCT" * 20 + "TAA"
    # only give sx and sy a CDS for OG_SC - sz's missing, so only 1 of the 3
    # possible pairs (sx vs sy) should be computable
    cds_by_protein = {"X_1": cds, "Y_1": cds}
    rows = run_pairwise_dnds(orthogroup_set, LOCUS_PREFIX_TO_SPECIES, cds_by_protein)
    og_sc_rows = [r for r in rows if r["orthogroup"] == "OG_SC"]
    assert len(og_sc_rows) == 1
    assert {og_sc_rows[0]["species_a"], og_sc_rows[0]["species_b"]} == {"sx", "sy"}


def test_run_pairwise_dnds_respects_max_orthogroups_cap():
    clusters = {f"OG{i}": [f"X_{i}", f"Y_{i}"] for i in range(10)}
    orthogroup_set = classify_orthogroups(clusters, {"X": "sx", "Y": "sy"}, ["sx", "sy"])
    cds = "ATG" + "GCT" * 10 + "TAA"
    cds_by_protein = {f"X_{i}": cds for i in range(10)}
    cds_by_protein.update({f"Y_{i}": cds for i in range(10)})
    rows = run_pairwise_dnds(orthogroup_set, {"X": "sx", "Y": "sy"}, cds_by_protein, max_orthogroups=3)
    assert len({r["orthogroup"] for r in rows}) == 3
