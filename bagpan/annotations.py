"""Parsers for bagRNA's current functional-annotation output.

bagRNA no longer runs `funannotate annotate`; its merge step
(bin/merge_functional_annotations.py, via modules/annotate_functional.nf) now
writes one gene-level `functional_annotation.tsv` covering everything
(product, GO, EC, KEGG, Pfam/InterPro, CAZy, MEROPS, PHI-base, secretion,
transmembrane, effector class, antiSMASH BGC role) - see that script's
`COLS` list, which parse_functional_annotation_tsv() mirrors exactly.

EffectorP3's own raw per-protein output is unchanged in format/location and
still carries a numeric probability the merged TSV drops (it only keeps the
class label), so parse_effectorp3() is kept as a supplementary source for
that one field.
"""

from __future__ import annotations

import dataclasses
import re
from pathlib import Path
from typing import Dict, Optional, Set, Tuple


@dataclasses.dataclass
class GeneFunctionalAnnotation:
    gene_id: str
    mrna_ids: Set[str]
    product: str
    gene_symbol: str
    go_terms: Set[str]
    ec_numbers: Set[str]
    kegg_ko: Set[str]
    kegg_pathways: Set[str]
    interpro: Set[str]
    pfam: Set[str]
    secreted: bool
    signalp: bool
    tm_tmbed: int
    tm_phobius: int
    sp_phobius: bool
    effector_class: str
    cazyme_family: str
    merops_hit: str
    merops_family: str
    rfam: Set[str]
    bgc_type: str
    bgc_role: str
    bgc_domains: Set[str]
    cog_category: str


def _split(value: str, sep: str) -> Set[str]:
    return {v.strip() for v in value.split(sep) if v.strip()}


_KEGG_PATHWAY_RE = re.compile(r"^(?:ko|map)(\d+)$")


def _split_kegg_pathways(value: str) -> Set[str]:
    """eggNOG's KEGG_Pathway field lists each pathway twice, once as a 'ko'
    (KEGG Orthology-based) id and once as the equivalent 'map' (reference
    pathway) id, e.g. 'ko00062|map00062' - both name the same pathway, so
    normalize both to a single canonical 'mapNNNNN' id.
    """
    ids: Set[str] = set()
    for token in value.split("|"):
        token = token.strip()
        if not token:
            continue
        match = _KEGG_PATHWAY_RE.match(token)
        ids.add(f"map{match.group(1)}" if match else token)
    return ids


def _to_int(value: str) -> int:
    return int(value) if value.strip().isdigit() else 0


def parse_functional_annotation_tsv(path: Optional[Path]) -> Dict[str, GeneFunctionalAnnotation]:
    """Parses functional_annotation.tsv into {gene_id: GeneFunctionalAnnotation}."""
    result: Dict[str, GeneFunctionalAnnotation] = {}
    if path is None:
        return result

    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}

        def get(row, name):
            i = col.get(name)
            return row[i] if i is not None and i < len(row) else ""

        for line in fh:
            if not line.strip():
                continue
            row = line.rstrip("\n").split("\t")
            gene_id = get(row, "gene_id")
            if not gene_id:
                continue
            result[gene_id] = GeneFunctionalAnnotation(
                gene_id=gene_id,
                mrna_ids=_split(get(row, "mrna_ids"), ","),
                product=get(row, "product"),
                gene_symbol=get(row, "gene_symbol"),
                go_terms=_split(get(row, "GO_terms"), "|"),
                ec_numbers=_split(get(row, "EC_numbers"), ","),
                kegg_ko=_split(get(row, "KEGG_KO"), ","),
                kegg_pathways=_split_kegg_pathways(get(row, "KEGG_pathways")),
                interpro=_split(get(row, "InterPro_accessions"), "|"),
                pfam=_split(get(row, "Pfam_domains"), "|"),
                secreted=get(row, "Secreted") == "Y",
                signalp=get(row, "SignalP") == "Y",
                tm_tmbed=_to_int(get(row, "TM_helices_TMbed")),
                tm_phobius=_to_int(get(row, "TM_helices_Phobius")),
                sp_phobius=get(row, "SP_Phobius") == "Y",
                effector_class=get(row, "EffectorP_class"),
                cazyme_family=get(row, "CAZyme_family"),
                merops_hit=get(row, "MEROPS_hit"),
                merops_family=get(row, "MEROPS_family"),
                rfam=_split(get(row, "Rfam_accessions"), "|"),
                bgc_type=get(row, "BGC_cluster_type"),
                bgc_role=get(row, "BGC_gene_role"),
                bgc_domains=_split(get(row, "BGC_domains"), "|"),
                cog_category=get(row, "COG_category"),
            )
    return result


def parse_annotation_stats(path: Optional[Path]) -> Dict[str, Tuple[int, str]]:
    """Parses bagRNA's per-species annotation_stats.txt
    ('Total genes                        17136  (100.0%)' style lines,
    written by bin/merge_functional_annotations.py) into
    {label: (count, pct_string)}.
    """
    result: Dict[str, Tuple[int, str]] = {}
    if path is None:
        return result
    pattern = re.compile(r"^(.+?)\s{2,}(\d+)\s+\(([\d.]+%)\)\s*$")
    with open(path) as fh:
        for line in fh:
            match = pattern.match(line.rstrip("\n"))
            if match:
                label, count, pct = match.groups()
                result[label.strip()] = (int(count), pct)
    return result


def parse_effectorp3(path: Optional[Path]) -> Dict[str, float]:
    """Parse EffectorP 3's raw output (unchanged by the funannotate removal):
    '# Identifier\\tCytoplasmic effector\\tApoplastic effector\\tNon-effector\\tPrediction'
    'SS109918_000002-T1 gene=... seq_id=... type=cds\\t-\\t-\\tY (0.984)\\tNon-effector'
    Returns {protein_id: probability} for anything not predicted
    'Non-effector'. Kept only for the numeric score - everything else about
    effector calls comes from functional_annotation.tsv's EffectorP_class.
    """
    effectors: Dict[str, float] = {}
    if path is None:
        return effectors
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            column = line.rstrip("\n").split("\t")
            if len(column) < 5:
                continue
            protein_id = column[0].split()[0]
            prediction = column[4].strip()
            if prediction == "Non-effector":
                continue
            score = None
            for field in column[1:4]:
                field = field.strip()
                if field.startswith("Y") and "(" in field:
                    score = float(field.split("(", 1)[1].rstrip(")"))
                    break
            effectors[protein_id] = score if score is not None else float("nan")
    return effectors
