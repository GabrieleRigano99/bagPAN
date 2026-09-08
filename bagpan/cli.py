from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List

from bagpan import dnds_runner, genes, go_enrichment, sequences, synteny
from bagpan.bigscape_runner import DEFAULT_IMAGE as DEFAULT_BIGSCAPE_IMAGE
from bagpan.bigscape_runner import BigscapeError, run_bigscape
from bagpan.categories import SpeciesAnnotations
from bagpan.discover import DiscoveryError, SpeciesFiles, discover_species
from bagpan.orthofinder_runner import DEFAULT_IMAGE as DEFAULT_ORTHOFINDER_IMAGE
from bagpan.orthofinder_runner import OrthoFinderError, run_orthofinder
from bagpan.orthogroups import OrthogroupsError, classify_orthogroups, parse_orthogroups_txt
from bagpan.report import run_report


def _species_arg(value: str) -> tuple:
    if "=" not in value:
        raise argparse.ArgumentTypeError(
            f"--species must be NAME=/path/to/bagRNA_outdir, got {value!r}"
        )
    name, path = value.split("=", 1)
    name, path = name.strip(), path.strip()
    if not name or not path:
        raise argparse.ArgumentTypeError(
            f"--species must be NAME=/path/to/bagRNA_outdir, got {value!r}"
        )
    return name, path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="bagpan",
        description=(
            "Comparative genomics / pangenome analysis over completed bagRNA "
            "output directories, using OrthoFinder orthogroups. Inspired by "
            "FunFinder_Pangenome.py (Wyka, 2019)."
        ),
    )
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="Run the full pangenome/comparative-genomics analysis")
    run.add_argument(
        "--species",
        action="append",
        type=_species_arg,
        required=True,
        metavar="NAME=PATH",
        help="A species/isolate: NAME=/path/to/completed/bagRNA/outdir. Repeatable, need >= 2.",
    )
    run.add_argument("--outdir", required=True, type=Path, help="Output directory")
    run.add_argument(
        "--skip-orthofinder",
        action="store_true",
        help="Don't run OrthoFinder; use --orthogroups instead",
    )
    run.add_argument(
        "--orthogroups",
        type=Path,
        help="Pre-computed Orthogroups.txt (required with --skip-orthofinder)",
    )
    run.add_argument(
        "--orthofinder-image",
        default=DEFAULT_ORTHOFINDER_IMAGE,
        help=f"Docker image:tag for OrthoFinder [default: {DEFAULT_ORTHOFINDER_IMAGE}]",
    )
    run.add_argument("--threads", type=int, default=4, help="Threads for OrthoFinder [default: 4]")
    run.add_argument(
        "--percent",
        type=float,
        default=0.5,
        help="Fraction of an orthogroup's species that must have a functional-category "
        "hit for the orthogroup to be classified as having that category "
        "(--classification-method threshold only) [default: 0.5]",
    )
    run.add_argument(
        "--classification-method",
        choices=["threshold", "mixture"],
        default="threshold",
        help="How to decide whether an orthogroup 'has' a functional category: a fixed "
        "--percent cutoff (default), or an opt-in 2-component Gaussian-mixture fit over "
        "each category's per-orthogroup hit fractions (PPanGGOLiN-flavored) [default: threshold]",
    )
    run.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="Alpha for BH-FDR correction of enrichment p-values [default: 0.05]",
    )
    run.add_argument(
        "--accumulation-permutations",
        type=int,
        default=200,
        help="Random species orderings averaged for the pangenome accumulation curve "
        "(exhaustive if <= n_species!) [default: 200]. Use 0 to skip the curve.",
    )
    run.add_argument(
        "--no-viz",
        action="store_true",
        help="Don't write the SVG plots / report.html (on by default)",
    )
    run.add_argument(
        "--run-bigscape",
        action="store_true",
        help="Also run BiG-SCAPE on each species' antiSMASH regions to cluster secondary-"
        "metabolite BGCs into gene-cluster families. Off by default: needs a large Pfam-A.hmm "
        "database and this integration hasn't been validated against a real run - see README.",
    )
    run.add_argument(
        "--bigscape-pfam",
        type=Path,
        help="Path to Pfam-A.hmm (required with --run-bigscape)",
    )
    run.add_argument(
        "--bigscape-image",
        default=DEFAULT_BIGSCAPE_IMAGE,
        help=f"Docker image:tag for BiG-SCAPE [default: {DEFAULT_BIGSCAPE_IMAGE}]",
    )
    run.add_argument(
        "--no-go-enrichment",
        action="store_true",
        help="Skip GO term enrichment (on by default)",
    )
    run.add_argument(
        "--go-obo",
        type=Path,
        help="Path to go-basic.obo. When given, GO annotations are propagated up the DAG "
        "(is_a/part_of) before testing, matching goatools/topGO. Without it, enrichment "
        "still runs but only on the exact annotated (leaf-level) terms.",
    )
    run.add_argument(
        "--go-min-count",
        type=int,
        default=3,
        help="Only test GO terms present in at least this many orthogroups pangenome-wide "
        "[default: 3]",
    )
    run.add_argument(
        "--run-dnds",
        action="store_true",
        help="Also compute pairwise dN/dS (Nei-Gojobori, pure Python - no mafft/PAML) on "
        "single-copy orthologs between every species pair. Off by default; requires "
        "--genome-fasta for every species.",
    )
    run.add_argument(
        "--genome-fasta",
        action="append",
        type=_species_arg,
        metavar="NAME=PATH",
        help="Genome assembly FASTA for a species (required per --species when --run-dnds "
        "is given) - used to splice CDS nucleotide sequences via annotated.gff3's coordinates.",
    )
    run.add_argument(
        "--dnds-max-orthogroups",
        type=int,
        default=200,
        help="Cap on the number of single-copy orthogroups analyzed for dN/dS [default: 200]",
    )

    return parser


def _run(args: argparse.Namespace) -> int:
    species_pairs: List[tuple] = args.species
    if len(species_pairs) < 2:
        print("ERROR: need at least 2 --species for a comparative analysis", file=sys.stderr)
        return 1

    names_seen = set()
    for name, _ in species_pairs:
        if name in names_seen:
            print(f"ERROR: duplicate --species name {name!r}", file=sys.stderr)
            return 1
        names_seen.add(name)

    if args.skip_orthofinder and not args.orthogroups:
        print("ERROR: --skip-orthofinder requires --orthogroups PATH", file=sys.stderr)
        return 1
    if args.run_bigscape and not args.bigscape_pfam:
        print("ERROR: --run-bigscape requires --bigscape-pfam PATH", file=sys.stderr)
        return 1

    genome_fasta_paths: Dict[str, str] = dict(args.genome_fasta or [])
    if args.run_dnds:
        missing = [name for name, _ in species_pairs if name not in genome_fasta_paths]
        if missing:
            print(
                f"ERROR: --run-dnds requires --genome-fasta NAME=PATH for every species; "
                f"missing: {', '.join(missing)}",
                file=sys.stderr,
            )
            return 1

    print(f"Discovering annotation files for {len(species_pairs)} species...")
    species_files: Dict[str, SpeciesFiles] = {}
    try:
        for name, path in species_pairs:
            sf = discover_species(name, path)
            species_files[name] = sf
            print(f"  {name}: proteome={sf.proteome} locus_prefix={sf.locus_prefix}")
    except DiscoveryError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    locus_prefix_to_species: Dict[str, str] = {}
    for name, sf in species_files.items():
        if sf.locus_prefix in locus_prefix_to_species:
            other = locus_prefix_to_species[sf.locus_prefix]
            print(
                f"ERROR: species {name!r} and {other!r} share the same locus-tag "
                f"prefix {sf.locus_prefix!r} - proteins won't be distinguishable",
                file=sys.stderr,
            )
            return 1
        locus_prefix_to_species[sf.locus_prefix] = name

    args.outdir.mkdir(parents=True, exist_ok=True)

    species_tree_path = statistics_path = None
    if args.skip_orthofinder:
        orthogroups_path = args.orthogroups
        orthofinder_cmd = None
    else:
        print(f"Running OrthoFinder ({args.orthofinder_image})...")
        try:
            of_result = run_orthofinder(
                list(species_files.values()), args.outdir, image=args.orthofinder_image, threads=args.threads
            )
        except OrthoFinderError as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            return 1
        orthogroups_path = of_result.orthogroups
        species_tree_path = of_result.species_tree
        statistics_path = of_result.statistics_overall
        orthofinder_cmd = f"docker run ... {args.orthofinder_image} orthofinder -f /data -t {args.threads}"
        print(f"  orthogroups: {orthogroups_path}")
        if species_tree_path:
            print(f"  species tree: {species_tree_path}")

    print("Parsing orthogroups...")
    try:
        clusters = parse_orthogroups_txt(orthogroups_path, locus_prefix_to_species)
        orthogroup_set = classify_orthogroups(clusters, locus_prefix_to_species, list(species_files.keys()))
    except OrthogroupsError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(f"  {len(clusters)} orthogroups")

    print("Parsing per-species functional annotations and gene structure...")
    species_annotations = {name: SpeciesAnnotations(sf) for name, sf in species_files.items()}
    gene_annotations = {name: genes.build_gene_annotation(sf) for name, sf in species_files.items()}

    protein_to_gene: Dict[str, str] = {}
    representative_transcripts = set()
    for ga in gene_annotations.values():
        protein_to_gene.update(ga.gene_of)
        representative_transcripts.update(ga.representative_of_gene.values())

    print("Checking orthogroup synteny support...")
    synteny_supported = synteny.check_synteny_support(orthogroup_set, gene_annotations)

    bigscape_result = None
    if args.run_bigscape:
        print(f"Running BiG-SCAPE ({args.bigscape_image})...")
        try:
            bigscape_result = run_bigscape(
                list(species_files.values()), args.outdir, args.bigscape_pfam,
                image=args.bigscape_image, threads=args.threads,
            )
        except BigscapeError as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            return 1
        print(f"  output: {bigscape_result.output_dir}")
        if not bigscape_result.gcf_by_bgc:
            print("  (no gene-cluster-family mapping could be resolved from BiG-SCAPE's output)")

    dnds_rows = None
    if args.run_dnds:
        print("Extracting CDS sequences and computing pairwise dN/dS...")
        cds_by_protein: Dict[str, str] = {}
        for name in species_files:
            genome = sequences.parse_fasta(genome_fasta_paths[name])
            cds_by_protein.update(dnds_runner.extract_all_cds(gene_annotations[name], genome))
        dnds_rows = dnds_runner.run_pairwise_dnds(
            orthogroup_set, locus_prefix_to_species, cds_by_protein,
            max_orthogroups=args.dnds_max_orthogroups,
        )
        n_orthogroups = len({row["orthogroup"] for row in dnds_rows})
        print(f"  {len(dnds_rows)} pairwise comparisons across {n_orthogroups} single-copy orthogroups")

    go_dag = None
    if not args.no_go_enrichment:
        if args.go_obo:
            print(f"Parsing GO DAG from {args.go_obo}...")
            go_dag = go_enrichment.parse_obo(args.go_obo)
        else:
            print("Running GO enrichment without term propagation (no --go-obo given)")

    print(f"Writing report to {args.outdir}...")
    run_meta = {
        "species": {
            name: {"bagrna_outdir": str(sf.root), "locus_prefix": sf.locus_prefix, "proteome": str(sf.proteome)}
            for name, sf in species_files.items()
        },
        "orthogroups_file": str(orthogroups_path),
        "orthofinder_command": orthofinder_cmd,
        "percent_threshold": args.percent,
        "classification_method": args.classification_method,
        "alpha": args.alpha,
        "dnds_max_orthogroups": args.dnds_max_orthogroups if args.run_dnds else None,
    }
    run_report(
        args.outdir,
        orthogroup_set,
        species_annotations,
        locus_prefix_to_species,
        protein_to_gene,
        representative_transcripts,
        synteny_supported,
        args.percent,
        args.alpha,
        args.classification_method,
        run_meta,
        species_tree_path=species_tree_path,
        statistics_path=statistics_path,
        accumulation_permutations=args.accumulation_permutations or None,
        make_viz=not args.no_viz,
        bigscape_result=bigscape_result,
        run_go_enrichment_flag=not args.no_go_enrichment,
        go_dag=go_dag,
        go_min_count=args.go_min_count,
        dnds_rows=dnds_rows,
    )
    print("Done.")
    return 0


def main(argv: List[str] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "run":
        return _run(args)
    parser.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())