"""Aggregate orthogroup classification + per-species functional annotations
into bagpan's output TSVs: pangenome stats, gene counts, per-category
orthogroup enrichment (Fisher's exact + BH-FDR, replicating the Core /
Accessory / Singleton comparison FunFinder_Pangenome.py's
write_final_results_and_figure() made, minus its plots), a combined
per-protein annotation table, the pangenome accumulation curve, and (when
available) OrthoFinder's own species tree/statistics and BiG-SCAPE's
gene-cluster-family assignments.
"""

from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Set

from bagpan import curves, go_enrichment, mixture, stats, viz
from bagpan.categories import CATEGORY_NAMES, SpeciesAnnotations
from bagpan.orthogroups import OrthogroupSet

PANGENOME_CLASSES = ("Core", "Accessory", "Singleton")


def _protein_species(protein_id: str, locus_prefix_to_species: Dict[str, str]) -> str:
    return locus_prefix_to_species[protein_id.split("_", 1)[0]]


def _global_category_sets(
    species_annotations: Dict[str, SpeciesAnnotations],
) -> Dict[str, Set[str]]:
    merged: Dict[str, Set[str]] = {name: set() for name in CATEGORY_NAMES}
    for sa in species_annotations.values():
        for name in CATEGORY_NAMES:
            merged[name] |= sa.categories[name]
    return merged


def write_pangenome_stats(outdir: Path, orthogroup_set: OrthogroupSet) -> Dict[str, int]:
    counts = {cls: 0 for cls in PANGENOME_CLASSES}
    for cluster in orthogroup_set.category.values():
        counts[cluster] += 1
    total = sum(counts.values())

    histogram: Dict[int, int] = {}
    for present in orthogroup_set.species_present.values():
        histogram[len(present)] = histogram.get(len(present), 0) + 1

    with open(outdir / "pangenome_stats.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["species_combinations", "n_orthogroups"])
        for n_species in sorted(histogram, reverse=True):
            w.writerow([n_species, histogram[n_species]])
        w.writerow([])
        for cls in PANGENOME_CLASSES:
            pct = (counts[cls] / total * 100) if total else 0.0
            w.writerow([cls, counts[cls], f"{pct:.2f}%"])
        w.writerow(["Total", total, "100.00%"])

    return counts


def write_gene_counts(
    outdir: Path,
    orthogroup_set: OrthogroupSet,
    locus_prefix_to_species: Dict[str, str],
    protein_to_gene: Dict[str, str],
) -> None:
    """Counts *distinct genes* per species per orthogroup (not raw
    transcript/protein rows) - a gene with several annotated isoforms must
    not inflate this count.
    """
    species_names = orthogroup_set.species_names
    with open(outdir / "gene_counts.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["orthogroup", *species_names])
        for cluster, proteins in orthogroup_set.proteins.items():
            genes_by_species: Dict[str, set] = {sp: set() for sp in species_names}
            for p in proteins:
                sp = _protein_species(p, locus_prefix_to_species)
                genes_by_species[sp].add(protein_to_gene.get(p, p))
            w.writerow([cluster, *(len(genes_by_species[sp]) for sp in species_names)])


def write_orthogroup_classification(
    outdir: Path,
    orthogroup_set: OrthogroupSet,
    protein_to_gene: Dict[str, str],
    locus_prefix_to_species: Dict[str, str],
    synteny_supported: Dict[str, Optional[bool]],
) -> None:
    with open(outdir / "orthogroup_classification.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            ["orthogroup", "pangenome_category", "n_species", "n_proteins", "n_genes", "synteny_supported", "species"]
        )
        for cluster, proteins in orthogroup_set.proteins.items():
            present = sorted(orthogroup_set.species_present[cluster])
            n_genes = len({protein_to_gene.get(p, p) for p in proteins})
            synteny_value = synteny_supported.get(cluster)
            synteny_str = "NA" if synteny_value is None else ("Yes" if synteny_value else "No")
            w.writerow(
                [cluster, orthogroup_set.category[cluster], len(present), len(proteins), n_genes, synteny_str, ",".join(present)]
            )


def write_functional_categories(
    outdir: Path,
    orthogroup_set: OrthogroupSet,
    species_annotations: Dict[str, SpeciesAnnotations],
    locus_prefix_to_species: Dict[str, str],
    protein_to_gene: Dict[str, str],
    percent: float,
    classification_method: str = "threshold",
) -> Dict[str, Dict[str, tuple]]:
    """Writes one TSV per functional category listing every orthogroup's hit
    status, and returns {category: {pangenome_class: (n_hit, n_total)}} for
    the enrichment summary.

    classification_method: "threshold" (default) calls an orthogroup a hit
    when >= `percent` of its species have the category; "mixture" instead
    fits a 2-component Gaussian mixture (bagpan.mixture) over the category's
    per-orthogroup hit-fractions and assigns hits by posterior probability -
    an opt-in, PPanGGOLiN-flavored statistical alternative to a fixed cutoff.
    """
    global_sets = _global_category_sets(species_annotations)
    cat_dir = outdir / "functional_category_orthogroups"
    cat_dir.mkdir(exist_ok=True)

    tallies: Dict[str, Dict[str, List[int]]] = {
        cat: {cls: [0, 0] for cls in PANGENOME_CLASSES} for cat in CATEGORY_NAMES
    }

    clusters = list(orthogroup_set.proteins)
    for category in CATEGORY_NAMES:
        category_set = global_sets[category]
        n_present_by_cluster = []
        hit_species_by_cluster = []
        fractions = []
        for cluster in clusters:
            proteins = orthogroup_set.proteins[cluster]
            present = orthogroup_set.species_present[cluster]
            n_present = len(present)
            hit_species = {
                _protein_species(p, locus_prefix_to_species)
                for p in proteins
                if protein_to_gene.get(p, p) in category_set
            }
            fraction = (len(hit_species) / n_present) if n_present else 0.0
            n_present_by_cluster.append(n_present)
            hit_species_by_cluster.append(hit_species)
            fractions.append(fraction)

        if classification_method == "mixture":
            hits = mixture.classify_by_mixture(fractions)
        else:
            hits = [f >= percent for f in fractions]

        with open(cat_dir / f"{category}.tsv", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(
                ["orthogroup", "pangenome_category", "n_species_in_cluster", "n_species_with_hit", "fraction_with_hit", "hit"]
            )
            for cluster, n_present, hit_species, fraction, is_hit in zip(
                clusters, n_present_by_cluster, hit_species_by_cluster, fractions, hits
            ):
                cls = orthogroup_set.category[cluster]
                tallies[category][cls][1] += 1
                if is_hit:
                    tallies[category][cls][0] += 1
                w.writerow([cluster, cls, n_present, len(hit_species), f"{fraction:.4f}", "Yes" if is_hit else "No"])

    return {cat: {cls: tuple(v) for cls, v in classes.items()} for cat, classes in tallies.items()}


def write_enrichment_summary(outdir: Path, tallies: Dict[str, Dict[str, tuple]], alpha: float) -> None:
    pairs = [("Core", "Accessory"), ("Core", "Singleton"), ("Accessory", "Singleton")]
    with open(outdir / "functional_enrichment_summary.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            [
                "category",
                "Core_hit", "Core_total", "Core_pct",
                "Accessory_hit", "Accessory_total", "Accessory_pct",
                "Singleton_hit", "Singleton_total", "Singleton_pct",
                "p_Core_vs_Accessory", "q_Core_vs_Accessory",
                "p_Core_vs_Singleton", "q_Core_vs_Singleton",
                "p_Accessory_vs_Singleton", "q_Accessory_vs_Singleton",
            ]
        )
        for category, classes in tallies.items():
            row = [category]
            for cls in PANGENOME_CLASSES:
                hit, total = classes[cls]
                pct = (hit / total * 100) if total else 0.0
                row += [hit, total, f"{pct:.2f}"]

            p_values = []
            for cls_a, cls_b in pairs:
                hit_a, total_a = classes[cls_a]
                hit_b, total_b = classes[cls_b]
                table = ((hit_a, total_a - hit_a), (hit_b, total_b - hit_b))
                p_values.append(stats.fisher_exact(table))
            q_values = stats.benjamini_hochberg(p_values)
            for p, q in zip(p_values, q_values):
                row += [f"{p:.6g}", f"{q:.6g}"]
            w.writerow(row)


def write_per_protein_annotations(
    outdir: Path,
    orthogroup_set: OrthogroupSet,
    species_annotations: Dict[str, SpeciesAnnotations],
    locus_prefix_to_species: Dict[str, str],
    protein_to_gene: Dict[str, str],
    representative_transcripts: Set[str],
) -> None:
    global_sets = _global_category_sets(species_annotations)
    with open(outdir / "per_protein_annotations.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            [
                "orthogroup", "cluster_size", "pangenome_category", "species", "gene_id",
                "protein_id", "is_representative_transcript",
                "gene_product", "pfam_domains", "iprscan_domains", "go_terms",
                "merops", "cazy", "secondary_metabolite_cluster",
                "secretome", "transmembrane", "effector", "effector_class", "effector_score",
            ]
        )
        for cluster, proteins in orthogroup_set.proteins.items():
            cluster_size = len(proteins)
            pan_category = orthogroup_set.category[cluster]
            for protein_id in proteins:
                species = _protein_species(protein_id, locus_prefix_to_species)
                sa = species_annotations[species]
                gene_id = protein_to_gene.get(protein_id, protein_id)
                is_secretome = gene_id in global_sets["secretome"]
                is_transmembrane = gene_id in global_sets["transmembrane"]
                is_effector = gene_id in global_sets["effectors"]
                is_metabolite = gene_id in global_sets["secondary_metabolites"]
                effector_score = sa.effector_scores.get(protein_id, "")
                w.writerow(
                    [
                        cluster, cluster_size, pan_category, species,
                        gene_id, protein_id,
                        "Yes" if protein_id in representative_transcripts else "No",
                        sa.product_of(gene_id), sa.pfam_of(gene_id),
                        sa.iprscan_of(gene_id), sa.go_of(gene_id),
                        sa.merops_of(gene_id), sa.cazy_of(gene_id),
                        "Yes" if is_metabolite else "No",
                        "Yes" if is_secretome else "No",
                        "Yes" if is_transmembrane else "No",
                        "Yes" if is_effector else "No",
                        sa.effector_class_of(gene_id),
                        effector_score,
                    ]
                )


def write_species_tree(outdir: Path, species_tree_path: Optional[Path]) -> None:
    if species_tree_path is not None and species_tree_path.is_file():
        shutil.copyfile(species_tree_path, outdir / "species_tree.nwk")


def _parse_orthofinder_statistics(path: Path) -> Dict[str, str]:
    stats_dict: Dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            if "\t" not in line:
                continue
            key, value = line.rstrip("\n").split("\t", 1)
            stats_dict[key] = value
    return stats_dict


def write_accumulation_curve(
    outdir: Path, orthogroup_set: OrthogroupSet, n_permutations: int
) -> tuple:
    curve = curves.accumulation_curve(orthogroup_set, n_permutations=n_permutations)
    ks = sorted(curve)

    with open(outdir / "pangenome_accumulation.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["genomes_sampled", "core_mean", "core_std", "pan_mean", "pan_std"])
        for k in ks:
            row = curve[k]
            w.writerow([k, f"{row['core_mean']:.3f}", f"{row['core_std']:.3f}", f"{row['pan_mean']:.3f}", f"{row['pan_std']:.3f}"])

    core_fit = curves.fit_power_law(ks, [curve[k]["core_mean"] for k in ks])
    pan_fit = curves.fit_power_law(ks, [curve[k]["pan_mean"] for k in ks])

    with open(outdir / "pangenome_openness.txt", "w") as fh:
        fh.write("Power-law fit y = kappa * k^gamma over the accumulation curve (k = genomes sampled).\n")
        fh.write("Descriptive diagnostic, not a hard open/closed classifier: a higher pan-genome\n")
        fh.write("gamma suggests a still-expanding pangenome; a gamma near 0 suggests a closed one.\n\n")
        if pan_fit:
            fh.write(f"pan-genome:  kappa={pan_fit[0]:.4f}  gamma={pan_fit[1]:.4f}\n")
        else:
            fh.write("pan-genome:  could not fit (insufficient data)\n")
        if core_fit:
            fh.write(f"core-genome: kappa={core_fit[0]:.4f}  gamma={core_fit[1]:.4f}\n")
        else:
            fh.write("core-genome: could not fit (insufficient data)\n")

    return curve, core_fit, pan_fit


def write_bigscape_output(outdir: Path, bigscape_result) -> None:
    if bigscape_result is None:
        return
    if bigscape_result.gcf_by_bgc:
        with open(outdir / "bigscape_gene_cluster_families.tsv", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["bgc_id", "gene_cluster_family"])
            for bgc_id, gcf in sorted(bigscape_result.gcf_by_bgc.items()):
                w.writerow([bgc_id, gcf])


def write_go_enrichment(
    outdir: Path,
    orthogroup_set: OrthogroupSet,
    species_annotations: Dict[str, SpeciesAnnotations],
    locus_prefix_to_species: Dict[str, str],
    protein_to_gene: Dict[str, str],
    percent: float,
    alpha: float,
    go_dag: Optional[go_enrichment.GoDag],
    go_min_count: int,
) -> None:
    # go_enrichment operates per protein id (it doesn't know about genes);
    # sa.go_terms is gene-keyed (functional_annotation.tsv is gene-level), so
    # resolve each orthogroup member's gene's GO set back onto its protein id.
    go_by_protein: Dict[str, Set[str]] = {}
    for cluster_proteins in orthogroup_set.proteins.values():
        for protein_id in cluster_proteins:
            species = _protein_species(protein_id, locus_prefix_to_species)
            gene_id = protein_to_gene.get(protein_id, protein_id)
            terms = species_annotations[species].go_terms.get(gene_id)
            if terms:
                go_by_protein[protein_id] = terms
    go_enrichment.run_go_enrichment(
        outdir, orthogroup_set, go_by_protein, locus_prefix_to_species, percent, alpha,
        dag=go_dag, min_population_count=go_min_count,
    )


def write_manifest(outdir: Path, meta: dict) -> None:
    with open(outdir / "run_manifest.json", "w") as fh:
        json.dump(meta, fh, indent=2, default=str)


def run_report(
    outdir: Path,
    orthogroup_set: OrthogroupSet,
    species_annotations: Dict[str, SpeciesAnnotations],
    locus_prefix_to_species: Dict[str, str],
    protein_to_gene: Dict[str, str],
    representative_transcripts: Set[str],
    synteny_supported: Dict[str, Optional[bool]],
    percent: float,
    alpha: float,
    classification_method: str,
    run_meta: dict,
    species_tree_path: Optional[Path] = None,
    statistics_path: Optional[Path] = None,
    accumulation_permutations: Optional[int] = 200,
    make_viz: bool = True,
    bigscape_result=None,
    run_go_enrichment_flag: bool = True,
    go_dag: Optional[go_enrichment.GoDag] = None,
    go_min_count: int = 3,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    pangenome_counts = write_pangenome_stats(outdir, orthogroup_set)
    write_gene_counts(outdir, orthogroup_set, locus_prefix_to_species, protein_to_gene)
    write_orthogroup_classification(
        outdir, orthogroup_set, protein_to_gene, locus_prefix_to_species, synteny_supported
    )
    tallies = write_functional_categories(
        outdir, orthogroup_set, species_annotations, locus_prefix_to_species, protein_to_gene,
        percent, classification_method,
    )
    write_enrichment_summary(outdir, tallies, alpha)
    write_per_protein_annotations(
        outdir, orthogroup_set, species_annotations, locus_prefix_to_species,
        protein_to_gene, representative_transcripts,
    )

    write_species_tree(outdir, species_tree_path)
    if statistics_path is not None and statistics_path.is_file():
        run_meta["orthofinder_statistics"] = _parse_orthofinder_statistics(statistics_path)

    curve = core_fit = pan_fit = None
    if accumulation_permutations:
        curve, core_fit, pan_fit = write_accumulation_curve(outdir, orthogroup_set, accumulation_permutations)

    write_bigscape_output(outdir, bigscape_result)

    if run_go_enrichment_flag:
        write_go_enrichment(
            outdir, orthogroup_set, species_annotations, locus_prefix_to_species, protein_to_gene,
            percent, alpha, go_dag, go_min_count,
        )
        run_meta["go_enrichment_propagated"] = go_dag is not None

    if make_viz:
        (outdir / "presence_absence_matrix.svg").write_text(viz.presence_absence_matrix_svg(orthogroup_set))
        (outdir / "category_enrichment.svg").write_text(viz.category_bar_svg(tallies))
        if curve:
            (outdir / "pangenome_accumulation.svg").write_text(
                viz.accumulation_curve_svg(curve, core_fit, pan_fit)
            )
        viz.write_report_html(outdir, orthogroup_set, pangenome_counts, tallies, curve, core_fit, pan_fit)

    write_manifest(outdir, run_meta)