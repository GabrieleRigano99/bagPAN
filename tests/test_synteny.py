from bagpan.genes import GeneAnnotation
from bagpan.orthogroups import classify_orthogroups
from bagpan.synteny import check_synteny_support

LOCUS_PREFIX_TO_SPECIES = {"X": "sx", "Y": "sy", "Z": "sz"}
SPECIES = ["sx", "sy", "sz"]

# OGA: syntenic - same flanking orthogroups (OG_L, OG_R) on both sides in both species
# OGB: non-syntenic - completely disjoint flanking orthogroups in each species
# OGC: singleton (sx only) -> nothing to compare -> NA
# OGD: sx + sz, but sz has no GFF3 gene-position data -> not comparable -> NA
CLUSTERS = {
    "OGA": ["X_1", "Y_1"],
    "OGB": ["X_2", "Y_2"],
    "OGC": ["X_3"],
    "OGD": ["X_4", "Z_1"],
    "OG_L": ["X_L", "Y_L"],
    "OG_R": ["X_R", "Y_R"],
    "OG_P": ["X_P"],
    "OG_Q": ["X_Q"],
    "OG_M": ["Y_M"],
    "OG_N": ["Y_N"],
}


def _gene_annotation(gene_positions: dict) -> GeneAnnotation:
    # gene_id == "transcript" id for this synthetic test: representative_of_gene is identity.
    return GeneAnnotation(
        gene_of={g: g for g in gene_positions},
        representative_of_gene={g: g for g in gene_positions},
        gene_position=gene_positions,
        from_gff3=True,
    )


def _build_scenario():
    orthogroup_set = classify_orthogroups(CLUSTERS, LOCUS_PREFIX_TO_SPECIES, SPECIES)

    sx = _gene_annotation({
        "X_L": ("c1", 100, "+"), "X_1": ("c1", 200, "+"), "X_R": ("c1", 300, "+"),
        "X_P": ("c2", 100, "+"), "X_2": ("c2", 200, "+"), "X_Q": ("c2", 300, "+"),
        "X_3": ("c3", 100, "+"),
        "X_4": ("c4", 100, "+"),
    })
    sy = _gene_annotation({
        "Y_L": ("c1", 100, "+"), "Y_1": ("c1", 200, "+"), "Y_R": ("c1", 300, "+"),
        "Y_M": ("c2", 100, "+"), "Y_2": ("c2", 200, "+"), "Y_N": ("c2", 300, "+"),
    })
    sz = _gene_annotation({})  # no GFF3 available for species z

    gene_annotations = {"sx": sx, "sy": sy, "sz": sz}
    return orthogroup_set, gene_annotations


def test_synteny_support_positive_case():
    orthogroup_set, gene_annotations = _build_scenario()
    result = check_synteny_support(orthogroup_set, gene_annotations, window=1)
    assert result["OGA"] is True


def test_synteny_support_negative_case():
    orthogroup_set, gene_annotations = _build_scenario()
    result = check_synteny_support(orthogroup_set, gene_annotations, window=1)
    assert result["OGB"] is False


def test_synteny_support_singleton_is_na():
    orthogroup_set, gene_annotations = _build_scenario()
    result = check_synteny_support(orthogroup_set, gene_annotations, window=1)
    assert result["OGC"] is None


def test_synteny_support_missing_gff3_species_is_na():
    orthogroup_set, gene_annotations = _build_scenario()
    result = check_synteny_support(orthogroup_set, gene_annotations, window=1)
    assert result["OGD"] is None