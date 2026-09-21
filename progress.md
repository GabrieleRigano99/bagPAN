# bagPAN — Progress Log

Comparative genomics / pangenome analysis tool for [bagRNA](https://github.com/GabrieleRigano99) annotation
outputs, using [OrthoFinder](https://github.com/davidemms/OrthoFinder) orthogroups.
Repo: https://github.com/GabrieleRigano99/bagPAN

Current state: **20 Python modules, ~3,400 lines, 123 passing tests, 9 commits**,
stdlib-only except OrthoFinder/BiG-SCAPE (Docker-shelled) and the optional
`--run-dnds` genome-FASTA input.

---

## 1. Initial build — core pangenome + functional enrichment (v1)

Inspired by `FunFinder_Pangenome.py` (Stephen A. Wyka, 2019 — the closest prior
art, built for funannotate+OrthoFinder output). Ported the core idea, dropped
everything requiring pandas/scipy/matplotlib/statsmodels/goatools:

- **`discover.py`** — locates each species' proteome + annotation files inside
  a bagRNA output dir; auto-detects the locus-tag prefix from FASTA headers
  rather than trusting a user-given name.
- **`orthogroups.py`** — parses OrthoFinder's `Orthogroups.txt`; classifies
  every orthogroup as **Core / Accessory / Singleton** using the species list
  explicitly given (not inferred from the largest cluster in the data, which
  is fragile).
- **`orthofinder_runner.py`** — stages proteomes, shells out to
  `docker run davidemms/orthofinder`.
- **`annotations.py` / `categories.py`** — parsed funannotate's per-tool
  output files (this whole layer was rewritten in step 4 below).
- **`stats.py`** — hand-rolled two-sided Fisher's exact test (hypergeometric
  enumeration via `math.comb`) and Benjamini-Hochberg FDR correction. No
  scipy/statsmodels. Cross-checked against hand-computed hypergeometric
  values and a known textbook example.
- **`report.py` / `cli.py`** — 7 functional categories (secretome,
  transmembrane, effectors, conserved domains, CAZymes, MEROPS peptidases,
  secondary metabolites) rolled up to orthogroup level via a `--percent`
  threshold, Fisher/BH-FDR tested Core vs Accessory vs Singleton.

**Validated**: full pytest suite, plus a live end-to-end run with real Docker
OrthoFinder against real bagRNA output (two pseudo-species built from
`newrun_Sschenckii_1099-18`, relabeled to a second locus prefix).

## 2. Landscape comparison → v2 roadmap

Researched the current pangenomics/comparative-genomics field (Panaroo,
PPanGGOLiN, Roary, Pangloss, GET_HOMOLOGUES, BiG-SCAPE) and ranked concrete
gaps. Implemented the agreed roadmap:

- **Isoform-inflation bug fix** (`genes.py`) — a gene with 3 annotated
  transcripts was being counted 3x in orthogroup/gene counts. Now collapses
  to one representative transcript (longest CDS) per gene, with a fallback
  heuristic (`-Tn` suffix stripping) when no GFF3 is available. This was a
  **real bug caught in bagPAN's own smoke test**, not a hypothetical.
- **Surfaced OrthoFinder's own outputs** instead of discarding them:
  `species_tree.nwk` (its rooted STAG tree) and `orthofinder_statistics` in
  the manifest (`Statistics_Overall.tsv`) — both free byproducts of the same
  Docker run.
- **`curves.py`** — pangenome accumulation (rarefaction) curve + Heaps'-law
  power-law openness fit, averaged over random species orderings
  (exhaustive when N! is small, sampled otherwise). Pure Python.
- **`viz.py`** — hand-rolled inline SVG (no matplotlib): presence/absence
  matrix, category enrichment bars, accumulation curve, assembled into a
  local `report.html`.
- **`synteny.py`** — lightweight QC signal: does the genomic neighborhood of
  an orthogroup's members actually corroborate the OrthoFinder clustering?
  Flag-only (never splits/merges), not a PanOCT/PPanGGOLiN-grade
  reconstruction.
- **`mixture.py`** — opt-in 2-component Gaussian-mixture classifier
  (`--classification-method mixture`) as a statistical alternative to the
  fixed `--percent` cutoff, in the spirit of PPanGGOLiN's partitioning.
  Pure Python EM, no scipy/sklearn.
- **`bigscape_runner.py`** — opt-in BiG-SCAPE integration (`--run-bigscape`)
  for secondary-metabolite gene-cluster-family clustering, Docker-shelled
  like OrthoFinder. SQLite output parsed defensively (schema introspected at
  runtime). **Flagged as unvalidated** — needs a large Pfam-A.hmm database
  and real antiSMASH output that wasn't practical to pull for testing here.

**Validated**: unit tests for every new piece, plus a repeat live Docker/
OrthoFinder run against real Sschenckii data exercising items 1-6 (BiG-SCAPE
excluded from live validation per agreed scope).

## 3. GO term enrichment (`go_enrichment.py`)

- Per-term Fisher's exact test, **Core/Accessory/Singleton vs. the rest of
  the pangenome**, BH-FDR corrected within each class — the same
  study-vs-population design FunFinder used via goatools, reimplemented on
  bagPAN's own stdlib Fisher/BH.
- Optional DAG propagation via `--go-obo /path/to/go-basic.obo` (hand-rolled
  OBO parser: `is_a`, `part_of`, `alt_id`, obsolete-term handling). Without
  it, only leaf-level annotated terms are tested (clearly logged).
- `--go-min-count` filters out terms too rare to ever be significant.

## 4. Adapted to bagRNA dropping `funannotate annotate`

bagRNA's functional-annotation stage changed: it now merges every tool's
output into one gene-level `functional_annotation/functional_annotation.tsv`
plus a final `functional_annotation/annotated.gff3`
(`bin/merge_functional_annotations.py`), and the proteome comes from AGAT
(`structural_annotation/final_proteins.faa`) instead of funannotate.

- Rewired `discover.py`, `annotations.py`, `categories.py` around the new
  three files. **Net simplification** — one well-structured TSV replaced
  ~10 fragile per-tool funannotate-format parsers, and it's already
  gene-level (no more antiSMASH-ID-cleanup hacks the old format needed).
- `report.py` updated everywhere it resolved protein→gene for category/GO
  lookups.
- EffectorP3's raw output (unchanged in format/location) kept as a
  supplementary source purely for the numeric probability score the merged
  TSV drops.
- **Validated**: all fixtures rebuilt to the new layout, plus a live run
  against real current-format `newrun_Sschenckii_1099-18` output.

## 5. funannotate-compare-inspired comparative breakdowns (`comparative_stats.py`)

Whole-genome (not orthology-dependent) class breakdowns, mirroring
`funannotate compare`'s CAZy/MEROPS/COG/secondary-metabolite/TF summary
tables — tabulated directly from columns `functional_annotation.tsv` already
carries:

- CAZyme family class (GH/GT/PL/CE/CBM/AA)
- MEROPS peptidase catalytic class
- COG functional category (including eggNOG's multi-letter fields, e.g. `"GM"`)
- Secondary-metabolite BGC type (hybrid clusters like `"NRPS/T1PKS"` split correctly)
- **Transcription-factor domain family** — 37 fungal TF-associated InterPro
  domains, taken verbatim from funannotate's own `tf_interpro.txt`
- **KEGG pathway** — a column bagRNA already produced but bagPAN hadn't
  parsed (only `KEGG_KO` was read before); normalizes eggNOG's duplicate
  `ko`/`map` id pairs to one canonical id
- Per-species `annotation_stats_summary.tsv`, sourced directly from bagRNA's
  own `annotation_stats.txt`

Each gets a TSV matrix + stacked-bar SVG (`viz.stacked_bar_svg`).

**Validated against real Sschenckii data** at every step: 33 TF families
detected (dominated by the Zn(2)-Cys(6) binuclear cluster domain — the
textbook fungal-specific TF expansion), 379 KEGG pathways (top hits:
Metabolic pathways, secondary-metabolite/antibiotic biosynthesis, ribosome —
exactly the standard fungal genome profile).

## 6. Pairwise dN/dS (`sequences.py`, `dnds.py`, `dnds_runner.py`)

The biggest lift of the funannotate-compare-inspired items — needed new
input data bagRNA doesn't publish (CDS nucleotide sequences).

- **`sequences.py`** — FASTA parsing + strand-aware, multi-exon CDS splicing
  from a genome assembly using GFF3 coordinates.
- **`genes.py`** extended to capture per-transcript CDS intervals/contig/
  strand (previously only a summed length was tracked).
- **`dnds.py`** — Needleman-Wunsch protein alignment (linear gap penalty, a
  deliberate simplicity tradeoff for already-orthologous sequences) →
  codon back-translation → Nei-Gojobori (1986) synonymous/nonsynonymous
  site and difference counting → Jukes-Cantor correction. Pure Python, no
  mafft/PAML.
- **`dnds_runner.py`** — selects single-copy-in-both-species orthogroup
  pairs, orchestrates the pipeline, caps runtime via
  `--dnds-max-orthogroups`.
- New required input for this feature only: `--genome-fasta NAME=PATH`
  (opt-in with `--run-dnds`), since bagRNA doesn't publish CDS sequences.

**Validated**: the NG86 algorithm was cross-checked during development
against biopython's independent `Bio.codonalign` implementation (15
randomized sequence pairs, matched to 4+ decimal places — biopython is not a
bagPAN dependency, just a one-off correctness oracle). CDS splicing and full
CLI wiring covered by unit/end-to-end tests with hand-derived synonymous/
nonsynonymous mutations. **Known gap, honestly documented**: never run
against a real multi-exon gene from a real genome FASTA — the matching
genome for the GFF3 used elsewhere in testing wasn't available locally
(searched, confirmed absent).

## 7. GitHub repository

Created `https://github.com/GabrieleRigano99/bagPAN` (private, MIT-licensed).
Every feature above has been committed and pushed incrementally with
descriptive commit messages.

## 8. `--species-dir`: simpler multi-species input

User feedback: repeating `--species NAME=PATH` was clunky. Added
`--species-dir DIR` (mutually exclusive with `--species`): every immediate
subdirectory becomes one species (subdirectory name = species name),
skipping hidden entries and non-directories. `--species` still works for
cases where directory names don't match desired species names.

## 9. Real-world validation runs

Located real, current-format completed bagRNA runs beyond the original
Sschenckii/synthetic-fixture testing (some in `/mnt/newvolume/bagRNA_test/`,
not the local repo):

- **Malassezia globosa vs. Candida orthopsilosis** (2 species): 6,860
  orthogroups, 37.1% Core / 62.9% Singleton (no Accessory possible with only
  2 species) — sensible for two evolutionarily distant fungi (Basidiomycete
  vs. Ascomycete).
- **+ Cryptococcus depauperatus** (3 species): 10,396 orthogroups, 23.0%
  Core / 12.0% Accessory / 65.1% Singleton. Functional enrichment showed
  secretome/effectors/CAZymes significantly enriched in Accessory vs Core
  (p ≈ 1e-12 to 2e-4) — the classic pattern of niche-adaptation genes
  skewing toward the variable genome. Flagged one caveat honestly:
  OrthoFinder's own species tree grouped the two Basidiomycetes wrongly
  (STAG topology unreliable with only 3 divergent species and limited
  single-copy orthologs — not a bagPAN bug).

## 10. Graphical report fixes

User feedback: "the graphical part is off, it's really bad." Diagnosed
concretely against the real 3-species run before fixing anything:

- **`presence_absence_matrix.svg`**: was one full-size cell per orthogroup,
  no real cap → **12,180px wide**, unusable. Fixed to derive cell width from
  a target plot width (~900px) regardless of orthogroup count → **1,080px**
  on the same real run.
- **`stacked_bar_svg`** (used for all comparative breakdowns): one legend
  entry + segment per distinct class → KEGG's ~380 pathways produced a
  **122KB** chart with an unreadable legend. Fixed to cap at the top ~20
  classes by count, folding the rest into "Other" → **7.4KB** on the same
  data.
- **GO enrichment had no chart at all** → added `go_enrichment_bar_svg`
  (top terms by q-value, red when below `--alpha`).
- **New, directly answering "which pathways do they share and which
  don't"**: `comparative_stats.category_membership_breakdown()` +
  `category_overlap_table()` → `comparative/kegg_pathway_overlap.tsv`/`.svg`.
  On the real 3-species run: 332 KEGG pathways shared by all 3, 18 shared by
  2, 20 unique to one.
- **`report.html` overhaul**: previously linked only 8 files and embedded 3
  charts. Now embeds every chart (GO per class, all 7 comparative
  breakdowns including the new KEGG overlap) under labeled sections, each
  wrapped in a scrollable container so an oversized chart never breaks the
  page layout, and links every output file across
  pangenome/comparative/go_enrichment/dnds.

**Validated**: re-ran the real 3-species comparison (reusing the
already-computed orthogroups, no need to redo the OrthoFinder Docker step)
and confirmed the exact before/after numbers above.

---

## Design principles that held throughout

- **stdlib-only wherever feasible.** Fisher's exact test, BH-FDR, GO DAG
  propagation, the mixture-model classifier, the accumulation-curve
  power-law fit, and the entire dN/dS pipeline (alignment + Nei-Gojobori
  counting) are all hand-rolled. The only external tools are OrthoFinder and
  (optionally) BiG-SCAPE, both Docker-shelled exactly the same way.
- **Never trust a claim without checking the actual file on disk.** Several
  real bugs and format-compatibility issues were caught this way (the
  isoform-inflation bug, the antiSMASH doubled-locus-prefix quirk, the
  funannotate-format removal, the graphical report width blowups) —
  consistently by inspecting real bagRNA output or real tool output before
  writing or fixing code, not by assumption.
- **Validate every non-trivial algorithm against an independent reference**
  when one exists (biopython for dN/dS) without adopting it as a runtime
  dependency.
- **Flag limitations honestly instead of overclaiming.** BiG-SCAPE's schema
  parsing and dN/dS-against-a-real-multi-exon-gene are both explicitly
  documented as unvalidated in the README, not silently assumed correct.

## Known gaps / possible future work

- Pfam/InterPro NMDS ordination (from `funannotate compare`) — not implemented.
- BiG-SCAPE integration unvalidated against a real run (needs a large
  Pfam-A.hmm database + real antiSMASH output).
- dN/dS pipeline unvalidated against a real multi-exon gene from a real
  genome FASTA (algorithm and CDS splicing are independently validated;
  that exact combination isn't).
- Genome fluidity and protein-length statistics from the original
  FunFinder script were never ported (deliberately, per original scoping).
