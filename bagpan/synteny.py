"""Lightweight synteny-based QC flag for orthogroups: does the genomic
neighborhood context of an orthogroup's members actually corroborate it, or
is it plausibly a spurious OrthoFinder clustering (e.g. driven by a
promiscuous shared domain)?

This is deliberately a flag-only heuristic signal, not a PanOCT/PPanGGOLiN-
grade synteny reconstruction: no orthogroup is ever split or merged here,
and it needs GFF3-derived gene positions (see genes.py) to say anything at
all - species without a GFF3 simply contribute no synteny evidence and
orthogroups only involving them come back "NA".
"""

from __future__ import annotations

import itertools
from collections import defaultdict
from typing import Dict, Optional

from bagpan.genes import GeneAnnotation
from bagpan.orthogroups import OrthogroupSet


def check_synteny_support(
    orthogroup_set: OrthogroupSet,
    gene_annotations: Dict[str, GeneAnnotation],
    window: int = 2,
) -> Dict[str, Optional[bool]]:
    """Returns {cluster_id: True/False/None} - None ("NA") for Singletons
    (nothing to compare against) and for orthogroups where no species pair
    has GFF3-derived positions for both members.
    """
    protein_to_cluster: Dict[str, str] = {}
    for cluster, proteins in orthogroup_set.proteins.items():
        for protein_id in proteins:
            protein_to_cluster[protein_id] = cluster

    contig_order: Dict[str, Dict[str, list]] = {}
    gene_cluster_of: Dict[str, Dict[str, Optional[str]]] = {}
    gene_position_index: Dict[str, Dict[str, tuple]] = {}

    for species, ga in gene_annotations.items():
        by_contig = defaultdict(list)
        cluster_of_gene: Dict[str, Optional[str]] = {}
        for gene_id, (contig, start, _strand) in ga.gene_position.items():
            by_contig[contig].append((start, gene_id))
            representative = ga.representative_of_gene.get(gene_id)
            cluster_of_gene[gene_id] = protein_to_cluster.get(representative) if representative else None
        for contig in by_contig:
            by_contig[contig].sort()
            by_contig[contig] = [gene_id for _start, gene_id in by_contig[contig]]

        position_index = {}
        for contig, ordered_genes in by_contig.items():
            for idx, gene_id in enumerate(ordered_genes):
                position_index[gene_id] = (contig, idx)

        contig_order[species] = dict(by_contig)
        gene_cluster_of[species] = cluster_of_gene
        gene_position_index[species] = position_index

    # cluster -> species -> the gene in that species belonging to the cluster
    cluster_gene_in_species: Dict[str, Dict[str, str]] = defaultdict(dict)
    for species, cluster_of_gene in gene_cluster_of.items():
        for gene_id, cluster in cluster_of_gene.items():
            if cluster is not None:
                cluster_gene_in_species[cluster][species] = gene_id

    def neighbor_clusters(species: str, gene_id: str) -> set:
        if gene_id not in gene_position_index.get(species, {}):
            return set()
        contig, idx = gene_position_index[species][gene_id]
        ordered_genes = contig_order[species][contig]
        neighbor_ids = ordered_genes[max(0, idx - window):idx] + ordered_genes[idx + 1:idx + 1 + window]
        return {
            gene_cluster_of[species][g]
            for g in neighbor_ids
            if gene_cluster_of[species].get(g) is not None
        }

    result: Dict[str, Optional[bool]] = {}
    for cluster, present_species in orthogroup_set.species_present.items():
        if len(present_species) < 2:
            result[cluster] = None
            continue

        comparable_pairs = 0
        supported_pairs = 0
        for sp_a, sp_b in itertools.combinations(sorted(present_species), 2):
            gene_a = cluster_gene_in_species.get(cluster, {}).get(sp_a)
            gene_b = cluster_gene_in_species.get(cluster, {}).get(sp_b)
            if gene_a is None or gene_b is None:
                continue
            comparable_pairs += 1
            if neighbor_clusters(sp_a, gene_a) & neighbor_clusters(sp_b, gene_b):
                supported_pairs += 1

        result[cluster] = None if comparable_pairs == 0 else (supported_pairs / comparable_pairs) >= 0.5

    return result