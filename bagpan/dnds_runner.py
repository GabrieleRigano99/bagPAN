"""Orchestrates pairwise dN/dS across single-copy orthologs: finds
orthogroup x species-pair combinations where each species contributes
exactly one protein (paralog-free), splices each protein's CDS nucleotide
sequence out of its species' genome FASTA using the GFF3-derived CDS
coordinates already parsed by bagpan.genes, and runs
bagpan.dnds.pairwise_dn_ds on each pair.
"""

from __future__ import annotations

import itertools
from typing import Dict, List

from bagpan import dnds, sequences
from bagpan.genes import GeneAnnotation
from bagpan.orthogroups import OrthogroupSet


def select_single_copy_pairs(
    orthogroup_set: OrthogroupSet, locus_prefix_to_species: Dict[str, str]
) -> Dict[str, Dict[str, str]]:
    """{cluster_id: {species: protein_id}} for orthogroups where every
    present species contributes exactly one protein - the only orthogroups a
    pairwise dN/dS comparison is well-defined on (no paralog ambiguity about
    which copy to compare).
    """
    result: Dict[str, Dict[str, str]] = {}
    for cluster, proteins in orthogroup_set.proteins.items():
        by_species: Dict[str, List[str]] = {}
        for protein_id in proteins:
            species = locus_prefix_to_species[protein_id.split("_", 1)[0]]
            by_species.setdefault(species, []).append(protein_id)
        if len(by_species) >= 2 and all(len(v) == 1 for v in by_species.values()):
            result[cluster] = {sp: v[0] for sp, v in by_species.items()}
    return result


def extract_all_cds(gene_annotation: GeneAnnotation, genome: Dict[str, str]) -> Dict[str, str]:
    """{transcript_id: cds_sequence} for every transcript this species'
    GeneAnnotation has CDS coordinates for (i.e. it came from a real GFF3,
    not the no-GFF3 fallback) and whose contig is present in `genome`.
    """
    cds_by_transcript: Dict[str, str] = {}
    for transcript_id, intervals in gene_annotation.transcript_cds_intervals.items():
        contig = gene_annotation.transcript_contig.get(transcript_id)
        strand = gene_annotation.transcript_strand.get(transcript_id)
        if not intervals or contig is None or contig not in genome:
            continue
        cds_by_transcript[transcript_id] = sequences.extract_cds_sequence(genome, contig, intervals, strand)
    return cds_by_transcript


def run_pairwise_dnds(
    orthogroup_set: OrthogroupSet,
    locus_prefix_to_species: Dict[str, str],
    cds_by_protein: Dict[str, str],
    max_orthogroups: int = 200,
) -> List[dict]:
    """Runs pairwise dN/dS for up to `max_orthogroups` single-copy
    orthogroups (chosen deterministically by sorted cluster id), across
    every pair of species present together in each. Skips any pair missing a
    CDS sequence (e.g. no GFF3 or no matching contig in the genome FASTA)
    rather than failing the whole run.
    """
    single_copy = select_single_copy_pairs(orthogroup_set, locus_prefix_to_species)
    rows: List[dict] = []
    for cluster in sorted(single_copy)[:max_orthogroups]:
        species_to_protein = single_copy[cluster]
        for sp_a, sp_b in itertools.combinations(sorted(species_to_protein), 2):
            protein_a, protein_b = species_to_protein[sp_a], species_to_protein[sp_b]
            cds_a = cds_by_protein.get(protein_a)
            cds_b = cds_by_protein.get(protein_b)
            if not cds_a or not cds_b:
                continue
            result = dnds.pairwise_dn_ds(cds_a, cds_b)
            rows.append(
                {
                    "orthogroup": cluster,
                    "species_a": sp_a, "species_b": sp_b,
                    "protein_a": protein_a, "protein_b": protein_b,
                    "n_codons_compared": result.n_codons_compared,
                    "dS": result.dS, "dN": result.dN, "omega": result.omega,
                }
            )
    return rows
