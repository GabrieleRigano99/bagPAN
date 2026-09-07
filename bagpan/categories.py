"""Combine bagRNA's gene-level functional_annotation.tsv into per-gene
functional category sets (secretome, transmembrane, effectors, conserved
domains, CAZymes, peptidases, secondary metabolites) for one species.

Everything here operates at gene granularity now (functional_annotation.tsv
is one row per gene, isoforms already collapsed) - callers resolve a
protein/transcript id to its gene id (bagpan.genes) before looking anything
up here.
"""

from __future__ import annotations

from typing import Dict, Set

from bagpan import annotations as ann
from bagpan.discover import SpeciesFiles

CATEGORY_NAMES = (
    "secretome",
    "transmembrane",
    "effectors",
    "conserved_domains",
    "cazymes",
    "peptidases",
    "secondary_metabolites",
)


class SpeciesAnnotations:
    """All per-gene annotation/classification results for one species."""

    def __init__(self, files: SpeciesFiles):
        self.files = files

        self.genes = ann.parse_functional_annotation_tsv(files.functional_annotation_tsv)
        # supplementary: EffectorP3's raw per-protein output still carries a
        # numeric probability the merged TSV doesn't (it only keeps the class
        # label) - protein-id keyed, unlike everything else here.
        self.effector_scores = ann.parse_effectorp3(files.effectorp3)

        self.go_terms: Dict[str, Set[str]] = {gid: g.go_terms for gid, g in self.genes.items() if g.go_terms}

        secretome = {gid for gid, g in self.genes.items() if g.secreted}
        transmembrane = {gid for gid, g in self.genes.items() if g.tm_tmbed > 0 or g.tm_phobius > 0}
        effectors = {gid for gid, g in self.genes.items() if g.effector_class}
        conserved = {gid for gid, g in self.genes.items() if g.interpro or g.pfam}
        cazymes = {gid for gid, g in self.genes.items() if g.cazyme_family}
        peptidases = {gid for gid, g in self.genes.items() if g.merops_hit}
        secondary_metabolites = {gid for gid, g in self.genes.items() if g.bgc_type or g.bgc_role}

        self.categories: Dict[str, Set[str]] = {
            "secretome": secretome,
            "transmembrane": transmembrane,
            "effectors": effectors,
            "conserved_domains": conserved,
            "cazymes": cazymes,
            "peptidases": peptidases,
            "secondary_metabolites": secondary_metabolites,
        }

    def product_of(self, gene_id: str) -> str:
        gene = self.genes.get(gene_id)
        return gene.product if gene else ""

    def pfam_of(self, gene_id: str) -> str:
        gene = self.genes.get(gene_id)
        return ",".join(sorted(gene.pfam)) if gene else ""

    def iprscan_of(self, gene_id: str) -> str:
        gene = self.genes.get(gene_id)
        return ",".join(sorted(gene.interpro)) if gene else ""

    def go_of(self, gene_id: str) -> str:
        return ",".join(sorted(self.go_terms.get(gene_id, [])))

    def merops_of(self, gene_id: str) -> str:
        gene = self.genes.get(gene_id)
        return gene.merops_hit if gene else ""

    def cazy_of(self, gene_id: str) -> str:
        gene = self.genes.get(gene_id)
        return gene.cazyme_family if gene else ""

    def effector_class_of(self, gene_id: str) -> str:
        gene = self.genes.get(gene_id)
        return gene.effector_class if gene else ""
