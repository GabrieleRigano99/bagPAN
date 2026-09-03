"""Parse OrthoFinder's Orthogroups.txt and classify orthogroups into
Core / Accessory / Singleton pangenome categories.

Ports the dictionary-building logic of FunFinder_Pangenome.py's
create_ortho_dictionaries(), but derives the "total number of genomes" from
the species explicitly passed to bagpan rather than inferring it from the
largest cluster observed in the data (which is fragile when no single
orthogroup happens to contain every species).
"""

from __future__ import annotations

import dataclasses
from pathlib import Path
from typing import Dict, List, Sequence, Set, Tuple


class OrthogroupsError(RuntimeError):
    pass


@dataclasses.dataclass
class OrthogroupSet:
    proteins: Dict[str, List[str]]
    species_present: Dict[str, Set[str]]
    category: Dict[str, str]
    species_names: Tuple[str, ...]


def parse_orthogroups_txt(
    path: "str | Path", locus_prefix_to_species: Dict[str, str]
) -> Dict[str, List[str]]:
    """Parse a classic OrthoFinder Orthogroups.txt file into
    {cluster_id: [protein_id, ...]}. Accepts the ':'-, tab-, and
    space-delimited layouts different orthology tools (OrthoFinder, OrthoMCL,
    SiLiX) have used.
    """
    clusters: Dict[str, List[str]] = {}
    with open(path) as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                continue
            if ":" in line:
                cluster, rest = line.split(":", 1)
            elif "\t" in line:
                cluster, rest = line.split("\t", 1)
            else:
                cluster, rest = line.split(" ", 1)
            proteins = rest.replace("\t", " ").split()
            if not proteins:
                continue
            unknown = {p for p in proteins if p.split("_", 1)[0] not in locus_prefix_to_species}
            if unknown:
                sample = ", ".join(sorted(unknown)[:3])
                raise OrthogroupsError(
                    f"Orthogroup {cluster} contains protein id(s) with a locus prefix "
                    f"bagpan doesn't recognize (e.g. {sample}) - not among the species "
                    "passed on the command line"
                )
            clusters[cluster] = proteins
    if not clusters:
        raise OrthogroupsError(f"No orthogroups parsed from {path}")
    return clusters


def classify_orthogroups(
    clusters: Dict[str, List[str]],
    locus_prefix_to_species: Dict[str, str],
    species_names: Sequence[str],
) -> OrthogroupSet:
    all_species = set(species_names)
    species_present: Dict[str, Set[str]] = {}
    category: Dict[str, str] = {}

    for cluster, proteins in clusters.items():
        present = {locus_prefix_to_species[p.split("_", 1)[0]] for p in proteins}
        species_present[cluster] = present
        if present == all_species:
            category[cluster] = "Core"
        elif len(present) == 1:
            category[cluster] = "Singleton"
        else:
            category[cluster] = "Accessory"

    return OrthogroupSet(
        proteins=clusters,
        species_present=species_present,
        category=category,
        species_names=tuple(species_names),
    )