"""Combine the raw annotation parsers into per-protein functional category
sets for one species (secretome, transmembrane, effectors, conserved
domains, CAZymes, peptidases, secondary metabolites). Ports the set-algebra
FunFinder_Pangenome.py's main_analyses_loop() used to derive secretome /
transmembrane / conserved from the raw per-tool hits.
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
    """All per-protein annotation/classification results for one species."""

    def __init__(self, files: SpeciesFiles):
        self.files = files

        self.gene_products = ann.parse_generic_annotations(files.genes_products)
        self.pfam = ann.parse_generic_annotations(files.pfam)
        self.iprscan = ann.parse_generic_annotations(files.iprscan)
        self.go_terms = ann.parse_go(files.iprscan)
        self.dbcan = ann.parse_generic_annotations(files.dbcan)
        self.merops = ann.parse_generic_annotations(files.merops)

        self.phobius_secreted, self.phobius_transmembrane = ann.parse_phobius(files.phobius)
        self.signalp_secreted = ann.parse_signalp(files.signalp)
        self.effectors = ann.parse_effectorp3(files.effectorp3)
        self.antismash_clusters, self.antismash_smcog = ann.parse_antismash(
            [files.antismash, files.antismash_clusters]
        )

        secretome = (self.phobius_secreted | self.signalp_secreted) - set(
            self.phobius_transmembrane
        )
        secondary_metabolites = self.antismash_clusters - self.antismash_smcog

        self.categories: Dict[str, Set[str]] = {
            "secretome": secretome,
            "transmembrane": set(self.phobius_transmembrane),
            "effectors": set(self.effectors),
            "conserved_domains": set(self.pfam) | set(self.iprscan),
            "cazymes": set(self.dbcan),
            "peptidases": set(self.merops),
            "secondary_metabolites": secondary_metabolites,
        }

    def product_of(self, protein_id: str) -> str:
        return ",".join(self.gene_products.get(protein_id, []))

    def pfam_of(self, protein_id: str) -> str:
        return ",".join(self.pfam.get(protein_id, []))

    def iprscan_of(self, protein_id: str) -> str:
        return ",".join(self.iprscan.get(protein_id, []))

    def go_of(self, protein_id: str) -> str:
        return ",".join(sorted(self.go_terms.get(protein_id, [])))

    def merops_of(self, protein_id: str) -> str:
        return ",".join(self.merops.get(protein_id, []))

    def cazy_of(self, protein_id: str) -> str:
        return ",".join(self.dbcan.get(protein_id, []))