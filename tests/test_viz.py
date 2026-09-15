import xml.etree.ElementTree as ET

from bagpan.orthogroups import classify_orthogroups
from bagpan.curves import accumulation_curve, fit_power_law
from bagpan.viz import (
    accumulation_curve_svg,
    category_bar_svg,
    go_enrichment_bar_svg,
    membership_breakdown_svg,
    presence_absence_matrix_svg,
    stacked_bar_svg,
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


def test_presence_absence_matrix_svg_width_stays_bounded_at_large_scale():
    # regression test: this used to render one full-size cell per orthogroup
    # with no real cap, producing a ~12000px-wide unusable image at real
    # pangenome scale (thousands of orthogroups).
    many_clusters = {f"OG{i}": [f"X_{i}", f"Y_{i}"] for i in range(6000)}
    orthogroup_set = classify_orthogroups(many_clusters, LOCUS_PREFIX_TO_SPECIES, ["sx", "sy"])
    svg = presence_absence_matrix_svg(orthogroup_set, target_plot_width=900)
    _assert_well_formed_svg(svg)
    root = ET.fromstring(svg)
    width = float(root.attrib["width"])
    assert width < 3000  # nowhere near the old ~12000px blowup


def test_presence_absence_matrix_svg_truncation_note_only_when_truncated():
    svg_small = presence_absence_matrix_svg(_orthogroup_set(), max_orthogroups=100)
    assert "showing first" not in svg_small

    many_clusters = {f"OG{i}": [f"X_{i}"] for i in range(50)}
    orthogroup_set = classify_orthogroups(many_clusters, {"X": "sx"}, ["sx"])
    svg_truncated = presence_absence_matrix_svg(orthogroup_set, max_orthogroups=10)
    assert "showing first 10 of 50 orthogroups" in svg_truncated


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


def test_stacked_bar_svg_is_well_formed():
    matrix = {"sx": {"GH": 3, "CBM": 1}, "sy": {"GH": 1, "AA": 2}}
    svg = stacked_bar_svg(matrix)
    _assert_well_formed_svg(svg)
    assert "sx" in svg and "sy" in svg
    assert "GH" in svg and "CBM" in svg and "AA" in svg


def test_stacked_bar_svg_handles_empty_matrix():
    _assert_well_formed_svg(stacked_bar_svg({}))
    _assert_well_formed_svg(stacked_bar_svg({"sx": {}}))


def test_stacked_bar_svg_caps_many_classes_into_other():
    # regression test: this used to render one legend entry + stacked segment
    # per distinct class with no cap - unusable for e.g. ~380 KEGG pathways.
    matrix = {"sx": {f"map{i:05d}": (100 - i) for i in range(50)}}
    svg = stacked_bar_svg(matrix, max_classes=10)
    _assert_well_formed_svg(svg)
    assert "Other" in svg
    # the 9 highest-count classes (100..92) should be kept individually
    assert "map00000" in svg
    assert "map00008" in svg
    # a low-count class should have been folded into "Other" instead
    assert "map00049" not in svg


def test_membership_breakdown_svg_well_formed_and_labels():
    svg = membership_breakdown_svg({3: 5, 2: 2, 1: 10}, n_species=3, title="test title")
    _assert_well_formed_svg(svg)
    assert "Shared by all" in svg
    assert "Unique (1 species)" in svg
    assert "Shared by 2 species" in svg
    assert "test title" in svg


def test_membership_breakdown_svg_handles_empty_breakdown():
    _assert_well_formed_svg(membership_breakdown_svg({}, n_species=2))


def test_go_enrichment_bar_svg_well_formed_and_respects_top_n():
    rows = [
        (f"GO:{i:07d}", 5, 10, 8, 20, 0.01, 0.01 * (i + 1))
        for i in range(20)
    ]
    svg = go_enrichment_bar_svg(rows, top_n=5)
    _assert_well_formed_svg(svg)
    assert svg.count("<rect") == 5
    assert "GO:0000000" in svg
    assert "GO:0000019" not in svg  # beyond top_n


def test_go_enrichment_bar_svg_uses_names_when_available():
    rows = [("GO:0000001", 3, 10, 5, 20, 0.001, 0.002)]
    svg = go_enrichment_bar_svg(rows, names={"GO:0000001": "some function"})
    assert "some function" in svg


def test_go_enrichment_bar_svg_handles_empty_rows():
    _assert_well_formed_svg(go_enrichment_bar_svg([]))


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


def test_write_report_html_embeds_go_enrichment_and_comparative_sections(tmp_path):
    orthogroup_set = _orthogroup_set()
    curve = accumulation_curve(orthogroup_set, n_permutations=50, seed=0)
    ks = sorted(curve)
    pan_fit = fit_power_law(ks, [curve[k]["pan_mean"] for k in ks])
    core_fit = fit_power_law(ks, [curve[k]["core_mean"] for k in ks])

    go_dir = tmp_path / "go_enrichment"
    go_dir.mkdir()
    (go_dir / "go_enrichment_core.svg").write_text('<svg xmlns="http://www.w3.org/2000/svg"><text>GO CORE MARKER</text></svg>')
    (go_dir / "go_enrichment_core.tsv").write_text("go_id\n")

    comp_dir = tmp_path / "comparative"
    comp_dir.mkdir()
    (comp_dir / "cog_category_counts.svg").write_text('<svg xmlns="http://www.w3.org/2000/svg"><text>COG MARKER</text></svg>')
    (comp_dir / "cog_category_counts.tsv").write_text("species\n")

    write_report_html(tmp_path, orthogroup_set, {"Core": 1, "Accessory": 1, "Singleton": 2}, {}, curve, core_fit, pan_fit)

    content = (tmp_path / "report.html").read_text()
    assert "GO CORE MARKER" in content
    assert "COG MARKER" in content
    assert "COG category" in content  # pretty-printed section title
    assert 'href="go_enrichment/go_enrichment_core.tsv"' in content
    assert 'href="comparative/cog_category_counts.tsv"' in content