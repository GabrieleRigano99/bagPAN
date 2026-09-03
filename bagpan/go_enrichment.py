"""Gene Ontology term enrichment across Core/Accessory/Singleton orthogroups:
a per-term Fisher's exact test (pangenome class vs. the rest of the
pangenome) with BH-FDR correction within each class - the same
study-set-vs-population design FunFinder_Pangenome.py used via goatools,
reimplemented on bagpan's own stdlib Fisher's exact / BH-FDR (bagpan.stats)
so no goatools dependency is needed.

Term propagation up the GO DAG (a child term's genes also count toward its
ancestors, the standard true-path rule) is optional: pass a GoDag built by
parse_obo() to enable it. Without one, only the exact annotated (leaf-level)
terms are tested - simpler, but not directly comparable to goatools/topGO
output, which always propagates.
"""

from __future__ import annotations

import csv
import dataclasses
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from bagpan import stats

PANGENOME_CLASSES = ("Core", "Accessory", "Singleton")


@dataclasses.dataclass
class GoDag:
    parents: Dict[str, Set[str]]
    names: Dict[str, str]
    alt_to_primary: Dict[str, str]


def parse_obo(path: "str | Path") -> GoDag:
    """Parses a go-basic.obo file into is_a/part_of parent relationships,
    term names, and alt_id -> primary-id aliases. Obsolete terms are
    dropped.
    """
    parents: Dict[str, Set[str]] = {}
    names: Dict[str, str] = {}
    alt_to_primary: Dict[str, str] = {}

    state = {"id": None, "parents": set(), "name": None, "alts": [], "obsolete": False}

    def flush() -> None:
        if state["id"] and not state["obsolete"]:
            parents[state["id"]] = set(state["parents"])
            if state["name"]:
                names[state["id"]] = state["name"]
            for alt in state["alts"]:
                alt_to_primary[alt] = state["id"]
        state.update(id=None, parents=set(), name=None, alts=[], obsolete=False)

    in_term = False
    with open(path) as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if line == "[Term]":
                flush()
                in_term = True
                continue
            if line.startswith("[") and line != "[Term]":
                flush()
                in_term = False
                continue
            if not in_term or not line:
                continue
            if line.startswith("id:"):
                state["id"] = line[len("id:"):].strip()
            elif line.startswith("name:"):
                state["name"] = line[len("name:"):].strip()
            elif line.startswith("alt_id:"):
                state["alts"].append(line[len("alt_id:"):].strip())
            elif line.startswith("is_a:"):
                target = line[len("is_a:"):].split("!")[0].strip()
                if target:
                    state["parents"].add(target)
            elif line.startswith("relationship: part_of"):
                target = line[len("relationship: part_of"):].split("!")[0].strip()
                if target:
                    state["parents"].add(target)
            elif line.startswith("is_obsolete:") and "true" in line.lower():
                state["obsolete"] = True
        flush()

    return GoDag(parents=parents, names=names, alt_to_primary=alt_to_primary)


def _resolve(go_id: str, dag: GoDag) -> str:
    return dag.alt_to_primary.get(go_id, go_id)


def _ancestors(go_id: str, dag: GoDag, cache: Dict[str, Set[str]]) -> Set[str]:
    go_id = _resolve(go_id, dag)
    if go_id in cache:
        return cache[go_id]
    result: Set[str] = set()
    cache[go_id] = result  # placeholder breaks cycles, mutated in place below
    for parent in dag.parents.get(go_id, ()):
        result.add(parent)
        result |= _ancestors(parent, dag, cache)
    return result


def propagate_go_terms(
    terms_by_protein: Dict[str, Set[str]], dag: Optional[GoDag]
) -> Dict[str, Set[str]]:
    """Expands each protein's GO term set to include all ancestors. A no-op
    (returns the input as-is) when dag is None.
    """
    if dag is None:
        return terms_by_protein
    cache: Dict[str, Set[str]] = {}
    propagated: Dict[str, Set[str]] = {}
    for protein_id, terms in terms_by_protein.items():
        expanded = set()
        for term in terms:
            resolved = _resolve(term, dag)
            expanded.add(resolved)
            expanded |= _ancestors(resolved, dag, cache)
        propagated[protein_id] = expanded
    return propagated


def orthogroup_go_associations(
    orthogroup_set,
    go_by_protein: Dict[str, Set[str]],
    locus_prefix_to_species: Dict[str, str],
    percent: float,
) -> Dict[str, Set[str]]:
    """{orthogroup: set of associated GO ids} - a term is associated with an
    orthogroup when >= `percent` of the orthogroup's present species have
    >= 1 member protein carrying it (mirrors bagpan's other functional
    categories).
    """
    associations: Dict[str, Set[str]] = {}
    for cluster, proteins in orthogroup_set.proteins.items():
        present = orthogroup_set.species_present[cluster]
        n_present = len(present)
        if n_present == 0:
            associations[cluster] = set()
            continue
        species_by_term: Dict[str, Set[str]] = {}
        for protein_id in proteins:
            species = locus_prefix_to_species[protein_id.split("_", 1)[0]]
            for term in go_by_protein.get(protein_id, ()):
                species_by_term.setdefault(term, set()).add(species)
        associations[cluster] = {
            term for term, species_set in species_by_term.items() if len(species_set) / n_present >= percent
        }
    return associations


EnrichmentRow = Tuple[str, int, int, int, int, float, float]  # go_id, study_hit, study_total, pop_hit, pop_total, p, q


def compute_go_enrichment(
    orthogroup_set,
    associations: Dict[str, Set[str]],
    min_population_count: int = 3,
) -> Dict[str, List[EnrichmentRow]]:
    """Pure computation (no I/O): returns {pangenome_class: [rows]}, each row
    a Fisher's-exact-tested, BH-FDR-corrected GO term, sorted by p-value.
    """
    population_count: Dict[str, int] = {}
    for terms in associations.values():
        for term in terms:
            population_count[term] = population_count.get(term, 0) + 1
    tested_terms = sorted(t for t, count in population_count.items() if count >= min_population_count)

    total_orthogroups = len(orthogroup_set.proteins)
    class_members: Dict[str, Set[str]] = {cls: set() for cls in PANGENOME_CLASSES}
    for cluster, cls in orthogroup_set.category.items():
        class_members[cls].add(cluster)

    term_clusters: Dict[str, Set[str]] = {term: set() for term in tested_terms}
    for cluster, terms in associations.items():
        for term in terms:
            if term in term_clusters:
                term_clusters[term].add(cluster)

    results: Dict[str, List[EnrichmentRow]] = {}
    for cls in PANGENOME_CLASSES:
        study_clusters = class_members[cls]
        study_total = len(study_clusters)
        rest_total = total_orthogroups - study_total

        p_values = []
        partial_rows = []
        for term in tested_terms:
            clusters_with_term = term_clusters[term]
            study_hit = len(clusters_with_term & study_clusters)
            pop_hit = len(clusters_with_term)
            rest_hit = pop_hit - study_hit
            table = ((study_hit, study_total - study_hit), (rest_hit, rest_total - rest_hit))
            p = stats.fisher_exact(table)
            p_values.append(p)
            partial_rows.append((term, study_hit, study_total, pop_hit, total_orthogroups))

        q_values = stats.benjamini_hochberg(p_values)
        rows: List[EnrichmentRow] = [
            (term, study_hit, study_total, pop_hit, pop_total, p, q)
            for (term, study_hit, study_total, pop_hit, pop_total), p, q in zip(partial_rows, p_values, q_values)
        ]
        rows.sort(key=lambda r: r[5])
        results[cls] = rows

    return results


def write_go_enrichment(
    outdir: Path, results: Dict[str, List[EnrichmentRow]], alpha: float, dag: Optional[GoDag] = None
) -> None:
    go_dir = outdir / "go_enrichment"
    go_dir.mkdir(exist_ok=True)
    for cls, rows in results.items():
        with open(go_dir / f"go_enrichment_{cls.lower()}.tsv", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(
                ["go_id", "name", "study_count", "study_total", "population_count", "population_total", "p_value", "q_value", "enriched"]
            )
            for go_id, study_hit, study_total, pop_hit, pop_total, p, q in rows:
                name = dag.names.get(go_id, "") if dag else ""
                w.writerow([go_id, name, study_hit, study_total, pop_hit, pop_total, f"{p:.6g}", f"{q:.6g}", "Yes" if q <= alpha else "No"])


def run_go_enrichment(
    outdir: Path,
    orthogroup_set,
    go_by_protein: Dict[str, Set[str]],
    locus_prefix_to_species: Dict[str, str],
    percent: float,
    alpha: float,
    dag: Optional[GoDag] = None,
    min_population_count: int = 3,
) -> Dict[str, List[EnrichmentRow]]:
    propagated = propagate_go_terms(go_by_protein, dag)
    associations = orthogroup_go_associations(orthogroup_set, propagated, locus_prefix_to_species, percent)
    results = compute_go_enrichment(orthogroup_set, associations, min_population_count)
    write_go_enrichment(outdir, results, alpha, dag)
    return results