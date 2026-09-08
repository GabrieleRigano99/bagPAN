# bagPAN

Comparative genomics / pangenome analysis over completed [bagRNA](../bagRNA)
output directories, using [OrthoFinder](https://github.com/davidemms/OrthoFinder)
orthogroups. Inspired by `FunFinder_Pangenome.py` (Stephen A. Wyka, 2019), and
built out further to close specific gaps against the wider pangenomics
landscape (Panaroo, PPanGGOLiN, Roary, Pangloss, GET_HOMOLOGUES, BiG-SCAPE).

Given N completed bagRNA runs (each with finished functional annotation),
bagPAN:

1. Locates each species' protein FASTA (`structural_annotation/final_proteins.faa`),
   gene-level `functional_annotation/functional_annotation.tsv`, and final
   `functional_annotation/annotated.gff3` (bagRNA's own merge step -
   `bin/merge_functional_annotations.py` - already aggregates every tool's
   output there; bagRNA no longer runs `funannotate annotate`), and
   auto-detects each species' locus-tag prefix from its own protein headers.
2. Runs OrthoFinder (via Docker) on the collected proteomes, or accepts a
   pre-computed `Orthogroups.txt` with `--skip-orthofinder`. When it runs
   OrthoFinder itself, it also surfaces OrthoFinder's own rooted species tree
   and comparative-genomics statistics instead of discarding them.
3. Collapses each gene to its representative transcript (longest CDS, from
   the GFF3) before counting - a gene with several annotated isoforms is one
   gene, not N. Classifies every orthogroup as Core / Accessory / Singleton.
4. Flags each orthogroup's synteny support: do the genomic neighbors of its
   members (via the GFF3) actually corroborate the clustering, or does it
   look like a spurious OrthoFinder call? A lightweight heuristic QC signal,
   not a re-clustering - see `bagpan/synteny.py`.
5. Classifies each gene into functional categories - secretome,
   transmembrane, predicted effectors (EffectorP3 class, from
   `functional_annotation.tsv`), conserved domains (Pfam/InterPro), CAZymes
   (dbCAN), MEROPS peptidases, secondary-metabolite cluster members
   (antiSMASH BGC role) - and rolls that up to orthogroup-level hits, via
   either a fixed `--percent` threshold (default) or an opt-in 2-component
   Gaussian-mixture fit (`--classification-method mixture`, in the spirit of
   PPanGGOLiN's statistical partitioning).
6. Runs a Fisher's exact test (Core vs Accessory vs Singleton, BH-FDR
   corrected) per category, and writes a combined per-protein annotation
   table.
7. Computes a pangenome accumulation (rarefaction) curve and a Heaps'-law
   power-law fit of pangenome openness, averaged over random species
   orderings.
8. Writes a static local `report.html` (presence/absence matrix, per-category
   enrichment bars, accumulation curve, all inline SVG - no matplotlib) plus
   the individual `.svg` files. On by default, `--no-viz` to skip.
9. Optionally (`--run-bigscape`, off by default) runs BiG-SCAPE on each
   species' antiSMASH regions to cluster secondary-metabolite BGCs into
   gene-cluster families - see the **BiG-SCAPE caveat** below.
10. Runs GO term enrichment per pangenome class (Core / Accessory / Singleton
    vs. the rest of the pangenome, Fisher's exact + BH-FDR - the same
    study-vs-population design FunFinder used via goatools). On by default;
    `--no-go-enrichment` to skip. Give `--go-obo /path/to/go-basic.obo` to
    propagate annotations up the DAG (is_a/part_of) before testing, matching
    goatools/topGO; without it, only the exact annotated (leaf-level) terms
    are tested.
11. Writes whole-genome (not orthology-dependent) comparative breakdowns
    inspired by `funannotate compare`'s CAZy/MEROPS/COG/TF/secondary-metabolite
    summary tables: per-species gene counts by CAZyme family class
    (GH/GT/PL/CE/CBM/AA), MEROPS peptidase class, COG functional category,
    antiSMASH BGC type, and transcription-factor family (via a curated list
    of 37 fungal TF-associated InterPro domains, taken verbatim from
    funannotate's own `tf_interpro.txt`) - each with a stacked-bar SVG - plus
    a genome/annotation-stats comparison table sourced directly from bagRNA's
    own per-species `annotation_stats.txt`.

Everything above is stdlib-only (Fisher's exact test, BH-FDR, the
accumulation-curve power-law fit, the mixture-model classifier, and the GO
DAG parser/propagation are all hand-rolled - no
pandas/scipy/matplotlib/statsmodels/goatools). Genome fluidity and
protein-length statistics from the original FunFinder script are still not
implemented.

## Usage

```bash
bagpan run \
    --species sporothrix_1099-18=/path/to/bagRNA/newrun_Sschenckii_1099-18 \
    --species sporothrix_ss02=/path/to/bagRNA/other_completed_run \
    --outdir results/ \
    --threads 8

# or, with orthogroups already computed elsewhere:
bagpan run \
    --species sp1=/path/to/run1 --species sp2=/path/to/run2 \
    --outdir results/ \
    --skip-orthofinder --orthogroups /path/to/OrthoFinder/Results_X/Orthogroups/Orthogroups.txt

# opt-in statistical partitioning instead of a fixed threshold:
bagpan run --species ... --outdir results/ --classification-method mixture

# also cluster secondary-metabolite BGCs into gene-cluster families:
bagpan run --species ... --outdir results/ --run-bigscape --bigscape-pfam /path/to/Pfam-A.hmm

# GO enrichment with full DAG propagation:
bagpan run --species ... --outdir results/ --go-obo /path/to/go-basic.obo
```

Each `--species` path must be a bagRNA output directory whose
`functional_annotation/` stage (`ANNOTATE_FUNCTIONAL`) has completed. Running
OrthoFinder requires Docker (pulls `davidemms/orthofinder`, override with
`--orthofinder-image`).

### Input layout

Per `--species` bagRNA output directory, bagPAN reads:

| File | Required | Purpose |
|---|---|---|
| `structural_annotation/final_proteins.faa` | yes | proteome fed to OrthoFinder; locus-tag prefix auto-detected from its headers |
| `functional_annotation/functional_annotation.tsv` | yes | gene-level product/GO/EC/KEGG/Pfam/InterPro/CAZy/MEROPS/PHI-base/secretion/TM/effector/BGC annotations (`bin/merge_functional_annotations.py`'s output) |
| `functional_annotation/annotated.gff3` | no (falls back to a `-Tn`-suffix heuristic) | gene/transcript/CDS structure, for representative-transcript selection and synteny |
| `functional_annotation/effectorp3_output.txt` | no | supplementary: the one field the merged TSV drops - EffectorP3's numeric probability |
| `functional_annotation/annotation_stats.txt` | no | feeds `comparative/annotation_stats_summary.tsv` |
| `functional_annotation/antismash_output/*.region*.gbk` | no | only used by `--run-bigscape` |

This matches bagRNA's current functional-annotation stage
(`modules/annotate_functional.nf`), which replaced `funannotate annotate`.

Key flags: `--percent` (threshold-method cutoff, default 0.5),
`--classification-method {threshold,mixture}`, `--alpha` (BH-FDR alpha),
`--accumulation-permutations` (default 200, `0` to skip the curve),
`--no-viz`, `--run-bigscape` / `--bigscape-pfam` / `--bigscape-image`,
`--no-go-enrichment`, `--go-obo`, `--go-min-count` (default 3).

## Output (`--outdir`)

- `pangenome_stats.tsv` - Core/Accessory/Singleton counts and species-count histogram
- `gene_counts.tsv` - **gene**-level (not transcript-level) count per species per orthogroup
- `orthogroup_classification.tsv` - orthogroup -> pangenome category, `n_genes`, `synteny_supported`, species present
- `functional_category_orthogroups/*.tsv` - one file per functional category, per-orthogroup hit status
- `functional_enrichment_summary.tsv` - per category: % of Core/Accessory/Singleton with a hit, pairwise Fisher p/BH-q values
- `per_protein_annotations.tsv` - one row per protein: orthogroup, pangenome category, `gene_id`, `is_representative_transcript`, and all functional annotations
- `pangenome_accumulation.tsv` / `pangenome_openness.txt` - the rarefaction curve and its power-law fit
- `go_enrichment/go_enrichment_{core,accessory,singleton}.tsv` - per-class GO term enrichment (Fisher's exact + BH-FDR)
- `species_tree.nwk` - OrthoFinder's rooted species tree (when it ran OrthoFinder itself)
- `report.html`, `presence_absence_matrix.svg`, `category_enrichment.svg`, `pangenome_accumulation.svg`
- `comparative/{cazyme_family_counts,merops_class_counts,cog_category_counts,secondary_metabolite_type_counts,transcription_factor_domain_counts}.tsv` (+ matching `.svg` stacked-bar charts) - per-species gene counts by class, whole-genome (not per-orthogroup)
- `comparative/annotation_stats_summary.tsv` - bagRNA's own per-species `annotation_stats.txt` numbers, aligned side by side
- `bigscape_gene_cluster_families.tsv` - BGC id -> gene-cluster-family, only when `--run-bigscape` resolved a mapping (see caveat below)
- `run_manifest.json` - inputs, species -> locus-prefix mapping, parameters used, and (when it ran) OrthoFinder's own overall statistics

## BiG-SCAPE caveat

`--run-bigscape` shells out to BiG-SCAPE 2.0's first-party Docker image
(`ghcr.io/medema-group/big-scape`) the same way bagPAN shells out to
OrthoFinder, and needs a Pfam-A.hmm database (`--bigscape-pfam`). **Unlike
the OrthoFinder integration, this has not been validated against a real
run** - doing so needs a large Pfam-A.hmm database and real antiSMASH output,
which wasn't practical to pull and run for validation here. Its SQLite output
is parsed defensively (table/column names are introspected at runtime rather
than assumed), and it degrades to a clear warning instead of crashing if the
schema doesn't match what this code expects - but treat
`bigscape_gene_cluster_families.tsv` as unverified until you've run it once
for real and checked it against BiG-SCAPE's own output directly. The mapping
is also only BGC-id -> gene-cluster-family; connecting individual proteins/
orthogroups to a specific BGC still means cross-referencing bagRNA's own
antiSMASH output by hand.

## Development

```bash
pip install -e .
pytest
```