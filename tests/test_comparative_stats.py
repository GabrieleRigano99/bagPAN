from bagpan import comparative_stats as cs
from bagpan.annotations import GeneFunctionalAnnotation


def _gene(gene_id, **overrides):
    defaults = dict(
        gene_id=gene_id, mrna_ids=set(), product="", gene_symbol="",
        go_terms=set(), ec_numbers=set(), kegg_ko=set(), interpro=set(), pfam=set(),
        secreted=False, signalp=False, tm_tmbed=0, tm_phobius=0, sp_phobius=False,
        effector_class="", cazyme_family="", merops_hit="", merops_family="",
        rfam=set(), bgc_type="", bgc_role="", bgc_domains=set(), cog_category="",
    )
    defaults.update(overrides)
    return GeneFunctionalAnnotation(**defaults)


class _FakeSpeciesAnnotations:
    def __init__(self, genes):
        self.genes = {g.gene_id: g for g in genes}


def test_cazy_class_extraction():
    assert cs.cazy_class("GH95_e26") == "GH"
    assert cs.cazy_class("CBM13_e436") == "CBM"
    assert cs.cazy_class("AA8_e49") == "AA"
    assert cs.cazy_class("") == ""


def test_merops_class_extraction():
    assert cs.merops_class("S09.951") == "S"
    assert cs.merops_class("T02.004") == "T"
    assert cs.merops_class("") == ""


def test_bgc_types_splits_hybrid_clusters():
    assert cs.bgc_types("NRPS-like/NRPS-like") == {"NRPS-like"}
    assert cs.bgc_types("NRPS/T1PKS") == {"NRPS", "T1PKS"}
    assert cs.bgc_types("") == set()


def test_cog_categories_splits_multiletter_field():
    assert cs.cog_categories("G") == {"G"}
    assert cs.cog_categories("GM") == {"G", "M"}
    assert cs.cog_categories("") == set()


def test_cazyme_family_counts_per_species():
    species_annotations = {
        "sx": _FakeSpeciesAnnotations([
            _gene("g1", cazyme_family="GH95_e26"),
            _gene("g2", cazyme_family="GH5_e436"),
            _gene("g3", cazyme_family="CBM13_e436"),
            _gene("g4"),  # no CAZyme hit
        ]),
        "sy": _FakeSpeciesAnnotations([_gene("g5", cazyme_family="AA8_e49")]),
    }
    counts = cs.cazyme_family_counts(species_annotations)
    assert counts["sx"] == {"GH": 2, "CBM": 1}
    assert counts["sy"] == {"AA": 1}


def test_merops_family_counts_per_species():
    species_annotations = {
        "sx": _FakeSpeciesAnnotations([
            _gene("g1", merops_family="S09.951"),
            _gene("g2", merops_family="S53.007"),
            _gene("g3", merops_family="T02.004"),
        ]),
    }
    counts = cs.merops_family_counts(species_annotations)
    assert counts["sx"] == {"S": 2, "T": 1}


def test_cog_category_counts_handles_multiletter_field():
    species_annotations = {
        "sx": _FakeSpeciesAnnotations([
            _gene("g1", cog_category="G"),
            _gene("g2", cog_category="GM"),
        ]),
    }
    counts = cs.cog_category_counts(species_annotations)
    # g1 contributes to G, g2 contributes to both G and M
    assert counts["sx"] == {"G": 2, "M": 1}


def test_secondary_metabolite_type_counts_splits_hybrids():
    species_annotations = {
        "sx": _FakeSpeciesAnnotations([
            _gene("g1", bgc_type="NRPS"),
            _gene("g2", bgc_type="NRPS/T1PKS"),
            _gene("g3"),  # no BGC
        ]),
    }
    counts = cs.secondary_metabolite_type_counts(species_annotations)
    assert counts["sx"] == {"NRPS": 2, "T1PKS": 1}


def test_all_class_dictionaries_use_uppercase_single_or_short_codes():
    assert set(cs.COG_DESCRIPTIONS) == set("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    assert set(cs.MEROPS_CLASS_DESCRIPTIONS) <= set("ACGMNPSTU")
    assert "GH" in cs.CAZY_CLASS_DESCRIPTIONS and "AA" in cs.CAZY_CLASS_DESCRIPTIONS
