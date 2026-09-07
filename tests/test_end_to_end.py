import csv
import json
from pathlib import Path

from bagpan.cli import main

FIXTURES = Path(__file__).parent / "fixtures"


def _read_tsv(path: Path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_run_end_to_end_with_precomputed_orthogroups(tmp_path):
    outdir = tmp_path / "out"
    argv = [
        "run",
        "--species", f"species_a={FIXTURES / 'species_a'}",
        "--species", f"species_b={FIXTURES / 'species_b'}",
        "--species", f"species_c={FIXTURES / 'species_c'}",
        "--outdir", str(outdir),
        "--skip-orthofinder",
        "--orthogroups", str(FIXTURES / "Orthogroups.txt"),
    ]
    rc = main(argv)
    assert rc == 0

    for name in [
        "pangenome_stats.tsv",
        "gene_counts.tsv",
        "orthogroup_classification.tsv",
        "functional_enrichment_summary.tsv",
        "per_protein_annotations.tsv",
        "run_manifest.json",
        "pangenome_accumulation.tsv",
        "pangenome_openness.txt",
        "report.html",
        "presence_absence_matrix.svg",
        "category_enrichment.svg",
        "pangenome_accumulation.svg",
    ]:
        assert (outdir / name).exists(), f"missing output {name}"
    # skip-orthofinder mode has no OrthoFinder run to source a species tree from
    assert not (outdir / "species_tree.nwk").exists()

    # orthogroup classification: 1 Core, 1 Accessory, 2 Singleton (see fixtures/Orthogroups.txt)
    classification = {row["orthogroup"]: row for row in _read_tsv(outdir / "orthogroup_classification.tsv")}
    assert classification["OG0000001"]["pangenome_category"] == "Core"
    assert classification["OG0000002"]["pangenome_category"] == "Accessory"
    assert classification["OG0000003"]["pangenome_category"] == "Singleton"
    assert classification["OG0000004"]["pangenome_category"] == "Singleton"

    # every fixture gene has exactly one transcript, so gene-level and protein-level
    # counts must agree here (the isoform-collapsing fix is exercised precisely by
    # test_genes.py's multi-isoform cases; this just checks the n_genes column is wired)
    for cluster in classification.values():
        assert cluster["n_genes"] == cluster["n_proteins"]

    # singletons have nothing to compare synteny against -> NA
    assert classification["OG0000003"]["synteny_supported"] == "NA"
    assert classification["OG0000004"]["synteny_supported"] == "NA"
    assert classification["OG0000001"]["synteny_supported"] in {"Yes", "No"}
    assert classification["OG0000002"]["synteny_supported"] in {"Yes", "No"}

    # gene-level counts per species per orthogroup (also 1:1 with proteins in this fixture)
    gene_counts = {row["orthogroup"]: row for row in _read_tsv(outdir / "gene_counts.tsv")}
    assert gene_counts["OG0000001"] == {"orthogroup": "OG0000001", "species_a": "1", "species_b": "1", "species_c": "1"}
    assert gene_counts["OG0000004"] == {"orthogroup": "OG0000004", "species_a": "0", "species_b": "0", "species_c": "1"}

    # per-protein table has exactly one row per protein across all 3 species (3+2+2 = 7),
    # each carrying its gene_id and representative-transcript flag
    per_protein = _read_tsv(outdir / "per_protein_annotations.tsv")
    assert len(per_protein) == 7
    for row in per_protein:
        assert row["gene_id"]
        assert row["is_representative_transcript"] == "Yes"  # single-isoform genes in this fixture

    # functional category tables all exist
    cat_dir = outdir / "functional_category_orthogroups"
    for category in [
        "secretome", "transmembrane", "effectors", "conserved_domains",
        "cazymes", "peptidases", "secondary_metabolites",
    ]:
        assert (cat_dir / f"{category}.tsv").exists()

    conserved = {row["orthogroup"]: row for row in _read_tsv(cat_dir / "conserved_domains.tsv")}
    assert conserved["OG0000001"]["hit"] == "Yes"  # all 3 species have a pfam/iprscan hit
    assert conserved["OG0000002"]["hit"] == "No"

    secretome = {row["orthogroup"]: row for row in _read_tsv(cat_dir / "secretome.tsv")}
    assert secretome["OG0000001"]["hit"] == "Yes"  # 2/3 species secreted, >= 0.5 threshold
    assert secretome["OG0000002"]["hit"] == "Yes"  # 2/2 species secreted

    effectors = {row["orthogroup"]: row for row in _read_tsv(cat_dir / "effectors.tsv")}
    assert effectors["OG0000004"]["hit"] == "Yes"
    assert effectors["OG0000001"]["hit"] == "No"

    transmembrane = {row["orthogroup"]: row for row in _read_tsv(cat_dir / "transmembrane.tsv")}
    assert transmembrane["OG0000004"]["hit"] == "Yes"

    cazymes = {row["orthogroup"]: row for row in _read_tsv(cat_dir / "cazymes.tsv")}
    assert cazymes["OG0000003"]["hit"] == "Yes"

    # secondary_metabolites: no fixture gene has BGC_cluster_type/BGC_gene_role set
    metabolites = {row["orthogroup"]: row for row in _read_tsv(cat_dir / "secondary_metabolites.tsv")}
    assert metabolites["OG0000001"]["hit"] == "No"

    manifest = json.loads((outdir / "run_manifest.json").read_text())
    assert set(manifest["species"]) == {"species_a", "species_b", "species_c"}
    assert manifest["species"]["species_a"]["locus_prefix"] == "SPA"
    assert manifest["go_enrichment_propagated"] is False  # no --go-obo given

    go_dir = outdir / "go_enrichment"
    for cls in ("core", "accessory", "singleton"):
        assert (go_dir / f"go_enrichment_{cls}.tsv").exists()


def test_run_with_no_go_enrichment_skips_it(tmp_path):
    outdir = tmp_path / "out"
    argv = [
        "run",
        "--species", f"species_a={FIXTURES / 'species_a'}",
        "--species", f"species_b={FIXTURES / 'species_b'}",
        "--species", f"species_c={FIXTURES / 'species_c'}",
        "--outdir", str(outdir),
        "--skip-orthofinder",
        "--orthogroups", str(FIXTURES / "Orthogroups.txt"),
        "--no-go-enrichment",
    ]
    assert main(argv) == 0
    assert not (outdir / "go_enrichment").exists()
    manifest = json.loads((outdir / "run_manifest.json").read_text())
    assert "go_enrichment_propagated" not in manifest


def test_run_with_go_obo_propagates_and_records_it_in_manifest(tmp_path):
    obo_path = tmp_path / "go-basic.obo"
    obo_path.write_text(
        "[Term]\nid: GO:0003824\nname: catalytic activity\n\n"
        "[Term]\nid: GO:0003700\nname: DNA-binding transcription factor activity\n"
        "is_a: GO:0003824 ! catalytic activity\n"
    )
    outdir = tmp_path / "out"
    argv = [
        "run",
        "--species", f"species_a={FIXTURES / 'species_a'}",
        "--species", f"species_b={FIXTURES / 'species_b'}",
        "--species", f"species_c={FIXTURES / 'species_c'}",
        "--outdir", str(outdir),
        "--skip-orthofinder",
        "--orthogroups", str(FIXTURES / "Orthogroups.txt"),
        "--go-obo", str(obo_path),
        "--go-min-count", "1",
        "--percent", "0.3",
    ]
    assert main(argv) == 0
    manifest = json.loads((outdir / "run_manifest.json").read_text())
    assert manifest["go_enrichment_propagated"] is True

    core = _read_tsv(outdir / "go_enrichment" / "go_enrichment_core.tsv")
    go_ids = {row["go_id"] for row in core}
    # the fixture only annotates the leaf GO:0003700, but propagation should
    # also surface its parent GO:0003824 with a resolved name
    assert "GO:0003700" in go_ids
    assert "GO:0003824" in go_ids
    parent_row = next(r for r in core if r["go_id"] == "GO:0003824")
    assert parent_row["name"] == "catalytic activity"


def test_run_with_no_viz_skips_plots(tmp_path):
    outdir = tmp_path / "out"
    argv = [
        "run",
        "--species", f"species_a={FIXTURES / 'species_a'}",
        "--species", f"species_b={FIXTURES / 'species_b'}",
        "--species", f"species_c={FIXTURES / 'species_c'}",
        "--outdir", str(outdir),
        "--skip-orthofinder",
        "--orthogroups", str(FIXTURES / "Orthogroups.txt"),
        "--no-viz",
    ]
    assert main(argv) == 0
    assert not (outdir / "report.html").exists()
    assert not (outdir / "presence_absence_matrix.svg").exists()
    assert (outdir / "orthogroup_classification.tsv").exists()


def test_run_with_mixture_classification_method(tmp_path):
    outdir = tmp_path / "out"
    argv = [
        "run",
        "--species", f"species_a={FIXTURES / 'species_a'}",
        "--species", f"species_b={FIXTURES / 'species_b'}",
        "--species", f"species_c={FIXTURES / 'species_c'}",
        "--outdir", str(outdir),
        "--skip-orthofinder",
        "--orthogroups", str(FIXTURES / "Orthogroups.txt"),
        "--classification-method", "mixture",
    ]
    assert main(argv) == 0
    manifest = json.loads((outdir / "run_manifest.json").read_text())
    assert manifest["classification_method"] == "mixture"
    cat_dir = outdir / "functional_category_orthogroups"
    assert (cat_dir / "secretome.tsv").exists()


def test_run_accumulation_permutations_zero_skips_curve(tmp_path):
    outdir = tmp_path / "out"
    argv = [
        "run",
        "--species", f"species_a={FIXTURES / 'species_a'}",
        "--species", f"species_b={FIXTURES / 'species_b'}",
        "--species", f"species_c={FIXTURES / 'species_c'}",
        "--outdir", str(outdir),
        "--skip-orthofinder",
        "--orthogroups", str(FIXTURES / "Orthogroups.txt"),
        "--accumulation-permutations", "0",
    ]
    assert main(argv) == 0
    assert not (outdir / "pangenome_accumulation.tsv").exists()
    assert not (outdir / "pangenome_openness.txt").exists()


def test_run_rejects_fewer_than_two_species(tmp_path):
    outdir = tmp_path / "out"
    argv = [
        "run",
        "--species", f"species_a={FIXTURES / 'species_a'}",
        "--outdir", str(outdir),
        "--skip-orthofinder",
        "--orthogroups", str(FIXTURES / "Orthogroups.txt"),
    ]
    rc = main(argv)
    assert rc == 1


def test_run_rejects_duplicate_locus_prefix(tmp_path):
    outdir = tmp_path / "out"
    argv = [
        "run",
        "--species", f"species_a={FIXTURES / 'species_a'}",
        "--species", f"species_a_dup={FIXTURES / 'species_a'}",
        "--outdir", str(outdir),
        "--skip-orthofinder",
        "--orthogroups", str(FIXTURES / "Orthogroups.txt"),
    ]
    rc = main(argv)
    assert rc == 1