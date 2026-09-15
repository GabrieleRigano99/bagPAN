"""Per-species comparative breakdowns inspired by funannotate compare's
CAZyme/MEROPS/COG/secondary-metabolite class summaries: tabulate gene counts
by *class* (not just the presence/absence bagpan.categories works with)
directly from functional_annotation.tsv's already-parsed columns, plus a
per-species stats comparison table sourced from bagRNA's own
annotation_stats.txt (written by bin/merge_functional_annotations.py).

Unlike bagpan.categories (which classifies at the orthogroup level), this
operates directly over each species' full gene set - a whole-genome view,
not an orthology-dependent one.
"""

from __future__ import annotations

import re
from collections import defaultdict
from typing import Callable, Dict, Iterable, List, Set, Tuple

from bagpan.annotations import GeneFunctionalAnnotation
from bagpan.categories import SpeciesAnnotations

COG_DESCRIPTIONS = {
    "A": "RNA processing and modification",
    "B": "Chromatin structure and dynamics",
    "C": "Energy production and conversion",
    "D": "Cell cycle control, cell division, chromosome partitioning",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "G": "Carbohydrate transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "J": "Translation, ribosomal structure and biogenesis",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "M": "Cell wall/membrane/envelope biogenesis",
    "N": "Cell motility",
    "O": "Posttranslational modification, protein turnover, chaperones",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport and catabolism",
    "R": "General function prediction only",
    "S": "Function unknown",
    "T": "Signal transduction mechanisms",
    "U": "Intracellular trafficking, secretion, and vesicular transport",
    "V": "Defense mechanisms",
    "W": "Extracellular structures",
    "X": "Mobilome: prophages, transposons",
    "Y": "Nuclear structure",
    "Z": "Cytoskeleton",
}

MEROPS_CLASS_DESCRIPTIONS = {
    "A": "Aspartic peptidase",
    "C": "Cysteine peptidase",
    "G": "Glutamic peptidase",
    "M": "Metallopeptidase",
    "N": "Asparagine peptide lyase",
    "P": "Mixed peptidase",
    "S": "Serine peptidase",
    "T": "Threonine peptidase",
    "U": "Peptidase of unknown catalytic type",
}

CAZY_CLASS_DESCRIPTIONS = {
    "GH": "Glycoside Hydrolase",
    "GT": "GlycosylTransferase",
    "PL": "Polysaccharide Lyase",
    "CE": "Carbohydrate Esterase",
    "CBM": "Carbohydrate-Binding Module",
    "AA": "Auxiliary Activity",
}

# Curated fungal transcription-factor InterPro domains, taken verbatim from
# funannotate's own funannotate/config/tf_interpro.txt (nextgenusfs/funannotate,
# BSD-licensed) - the same list funannotate compare's transcription-factor
# comparison uses.
TF_INTERPRO_DESCRIPTIONS = {
    "IPR000967": "NF-X1-type zinc finger",
    "IPR006856": "Mating-type protein MAT alpha-1 HMG-Box",
    "IPR018501": "DDT domain",
    "IPR007196": "CCR4-Not complex component",
    "IPR007396": "Putative FMN-binding domain",
    "IPR004595": "TFIIH C1-like domain",
    "IPR004181": "MIZ zinc finger",
    "IPR000818": "TEA/ATTS domain family",
    "IPR001387": "Helix-turn-helix",
    "IPR001289": "CCAAT-binding TF (CBF-B/NF-YA) subunit B",
    "IPR003120": "STE-like TF",
    "IPR003150": "RFX DNA-binding domain",
    "IPR004198": "Zinc finger",
    "IPR007604": "CP2 TF",
    "IPR008895": "YL1 nuclear protein",
    "IPR010770": "SGT1 protein",
    "IPR018004": "KilA-N domain",
    "IPR024061": "NDT80/PhoG-like DNA-binding family",
    "IPR003656": "BED zinc finger",
    "IPR018060": "Bacterial regulatory HTH proteins",
    "IPR005011": "SART-1 family",
    "IPR002100": "SRF-type TF (DNA-binding and dimerization domain)",
    "IPR010666": "GRF zinc finger",
    "IPR000232": "HSF-type DNA-binding",
    "IPR001766": "Fork head domain",
    "IPR001356": "Homeobox domain",
    "IPR000679": "GATA zinc finger",
    "IPR013767": "PAS fold",
    "IPR001878": "Zinc knuckle (CCHC)",
    "IPR004827": "Basic region leucine zipper",
    "IPR007889": "Helix-turn-helix",
    "IPR011598": "Helix-loop-helix DNA-binding domain",
    "IPR001005": "Myb-like DNA-binding domain",
    "IPR007087": "Zinc finger (C2H2)",
    "IPR001138": "Fungal Zn(2)-Cys(6) binuclear cluster domain",
    "IPR007219": "Fungal-specific TF domain",
}

_CAZY_CLASS_RE = re.compile(r"^([A-Za-z]+)\d")


def cazy_class(family: str) -> str:
    """'GH95_e26' -> 'GH', 'CBM13_e436' -> 'CBM'."""
    match = _CAZY_CLASS_RE.match(family)
    return match.group(1).upper() if match else family


def merops_class(family: str) -> str:
    """'S09.951' -> 'S' (MEROPS catalytic-type letter)."""
    return family[0].upper() if family else ""


def bgc_types(bgc_type_field: str) -> Set[str]:
    """'NRPS-like/NRPS-like' (hybrid clusters, slash-joined) -> {'NRPS-like'}."""
    return {t.strip() for t in bgc_type_field.split("/") if t.strip()}


def cog_categories(cog_field: str) -> Set[str]:
    """eggNOG's COG_category can carry multiple letters in one field ('GM')."""
    return {c for c in cog_field.strip() if c.isalpha()}


def _count_by_class(
    species_annotations: Dict[str, SpeciesAnnotations],
    extractor: Callable[[GeneFunctionalAnnotation], Iterable[str]],
) -> Dict[str, Dict[str, int]]:
    """extractor(gene) -> class labels this gene contributes to (usually 0 or
    1, sometimes >1 for hybrid BGCs / multi-letter COG fields).
    """
    counts: Dict[str, Dict[str, int]] = {}
    for species, sa in species_annotations.items():
        tally: Dict[str, int] = {}
        for gene in sa.genes.values():
            for label in extractor(gene):
                tally[label] = tally.get(label, 0) + 1
        counts[species] = tally
    return counts


def cazyme_family_counts(species_annotations: Dict[str, SpeciesAnnotations]) -> Dict[str, Dict[str, int]]:
    return _count_by_class(
        species_annotations, lambda g: [cazy_class(g.cazyme_family)] if g.cazyme_family else []
    )


def merops_family_counts(species_annotations: Dict[str, SpeciesAnnotations]) -> Dict[str, Dict[str, int]]:
    return _count_by_class(
        species_annotations, lambda g: [merops_class(g.merops_family)] if g.merops_family else []
    )


def cog_category_counts(species_annotations: Dict[str, SpeciesAnnotations]) -> Dict[str, Dict[str, int]]:
    return _count_by_class(species_annotations, lambda g: cog_categories(g.cog_category))


def secondary_metabolite_type_counts(species_annotations: Dict[str, SpeciesAnnotations]) -> Dict[str, Dict[str, int]]:
    return _count_by_class(species_annotations, lambda g: bgc_types(g.bgc_type))


def transcription_factor_domains(gene: GeneFunctionalAnnotation) -> Set[str]:
    return gene.interpro & TF_INTERPRO_DESCRIPTIONS.keys()


def transcription_factor_domain_counts(species_annotations: Dict[str, SpeciesAnnotations]) -> Dict[str, Dict[str, int]]:
    return _count_by_class(species_annotations, transcription_factor_domains)


def kegg_pathway_counts(species_annotations: Dict[str, SpeciesAnnotations]) -> Dict[str, Dict[str, int]]:
    """Per-species gene counts by KEGG pathway (canonical 'mapNNNNN' ids -
    see annotations._split_kegg_pathways for the ko/map dedup)."""
    return _count_by_class(species_annotations, lambda g: g.kegg_pathways)


def _label_species_sets(matrix: Dict[str, Dict[str, int]]) -> Dict[str, Set[str]]:
    label_species: Dict[str, Set[str]] = defaultdict(set)
    for species, tally in matrix.items():
        for label, count in tally.items():
            if count > 0:
                label_species[label].add(species)
    return label_species


def category_membership_breakdown(matrix: Dict[str, Dict[str, int]]) -> Dict[int, int]:
    """{n_species_sharing: n_labels} - e.g. {3: 120, 2: 45, 1: 800} means 120
    labels (KEGG pathways, CAZy families, ...) are present in all 3 species,
    45 in exactly 2, and 800 are unique to a single species. Answers "how
    much is shared vs. species-specific" directly, from any of the
    comparative_stats *_counts() matrices.
    """
    breakdown: Dict[int, int] = {}
    for species in _label_species_sets(matrix).values():
        breakdown[len(species)] = breakdown.get(len(species), 0) + 1
    return breakdown


def category_overlap_table(matrix: Dict[str, Dict[str, int]]) -> List[Tuple[str, int, List[str]]]:
    """[(label, n_species, [species, ...])], sorted by n_species descending
    then label - the per-label detail behind category_membership_breakdown().
    """
    rows = [(label, len(species), sorted(species)) for label, species in _label_species_sets(matrix).items()]
    rows.sort(key=lambda row: (-row[1], row[0]))
    return rows
