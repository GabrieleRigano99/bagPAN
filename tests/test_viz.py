import xml.etree.ElementTree as ET

from bagpan.orthogroups import classify_orthogroups
from bagpan.curves import accumulation_curve, fit_power_law
from bagpan.viz import (
    accumulation_curve_svg,
    category_bar_svg,
    presence_absence_matrix_svg,
    write_report_html,
)

LOCUS_PREFIX_TO_SPECIES = {"X": "sx", "Y": "sy", "Z": "sz"}
SPECIES = ["sx", "sy", "sz"]
CLUSTERS = {
    "OG1": ["X_1", "Y_1", "Z_1"],
    "OG2": ["X_2", "Y_2"],
    "OG3": ["X_3"],
    "OG4": ["Z_2"],
}


def _orthogroup_set():
    return classify_orthogroups(CLUSTERS, LOCUS_PREFIX_TO_SPECIES, SPECIES)


def _assert_well_formed_svg(svg_text: str) -> None:
    root = ET.fromstring(svg_text)  # raises ParseError on malformed XML
    assert root.tag.endswith("svg")


def test_presence_absence_matrix_svg_is_well_formed_and_has_species_labels():
    svg = presence_absence_matrix_svg(_orthogroup_set())
    _assert_well_formed_svg(svg)
    for sp in SPECIES:
        assert sp in svg


def test_category_bar_svg_is_well_formed():
    tallies = {
        "secretome": {"Core": (1, 1), "Accessory": (1, 2), "Singleton": (0, 2)},
        "cazymes": {"Core": (0, 1), "Accessory": (0, 2), "Singleton": (1, 2)},
    }
    svg = category_bar_svg(tallies)
    _assert_well_formed_svg(svg)
    assert "secretome" in svg and "cazymes" in svg


def test_accumulation_curve_svg_is_well_formed_with_and_without_fit():
    curve = accumulation_curve(_orthogroup_set(), n_permutations=50, seed=0)
    svg_no_fit = accumulation_curve_svg(curve)
    _assert_well_formed_svg(svg_no_fit)

    ks = sorted(curve)
    pan_fit = fit_power_law(ks, [curve[k]["pan_mean"] for k in ks])
    core_fit = fit_power_law(ks, [curve[k]["core_mean"] for k in ks])
    svg_with_fit = accumulation_curve_svg(curve, core_fit=core_fit, pan_fit=pan_fit)
    _assert_well_formed_svg(svg_with_fit)


def test_accumulation_curve_svg_handles_empty_curve():
    _assert_well_formed_svg(accumulation_curve_svg({}))


def test_write_report_html(tmp_path):
    orthogroup_set = _orthogroup_set()
    (tmp_path / "pangenome_stats.tsv").write_text("dummy\n")
    tallies = {"secretome": {"Core": (1, 1), "Accessory": (0, 1), "Singleton": (0, 2)}}
    curve = accumulation_curve(orthogroup_set, n_permutations=50, seed=0)
    ks = sorted(curve)
    pan_fit = fit_power_law(ks, [curve[k]["pan_mean"] for k in ks])
    core_fit = fit_power_law(ks, [curve[k]["core_mean"] for k in ks])

    write_report_html(tmp_path, orthogroup_set, {"Core": 1, "Accessory": 1, "Singleton": 2}, tallies, curve, core_fit, pan_fit)

    report = tmp_path / "report.html"
    assert report.exists()
    content = report.read_text()
    assert "<svg" in content
    assert 'href="pangenome_stats.tsv"' in content
    # a TSV that wasn't created shouldn't get a dangling link
    assert 'href="gene_counts.tsv"' not in content