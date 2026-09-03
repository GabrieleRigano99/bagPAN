"""Per-format parsers for the annotation files bagRNA (via funannotate)
produces. Where the file format matches funannotate's native output
(pfam/iprscan/dbCAN/merops/genes-products/antismash/phobius/signalp), these
port FunFinder_Pangenome.py's parsers directly. EffectorP3's column layout
differs from the EffectorP2 format FunFinder was written against, so it gets
a new parser.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple


def parse_generic_annotations(path: Optional[Path]) -> Dict[str, List[str]]:
    """Parse a funannotate-style 'annotations.<db>.txt' file:
    'protein_id\\tfield\\tvalue' per line, e.g.
        SS109918_014817-T2   db_xref   PFAM:PF00198
        SS109918_013908-T1   product   Putative aryl-alcohol dehydrogenase aad14
    GO term lines (field containing 'go') are skipped here - use parse_go().
    'name' lines are skipped. For db_xref/note-style lines the value after
    the first ':' is kept (e.g. 'PF00198' from 'PFAM:PF00198'); lines with no
    ':' are kept as-is (e.g. free-text products).
    """
    result: Dict[str, List[str]] = {}
    if path is None:
        return result
    with open(path) as fh:
        for line in fh:
            column = line.rstrip("\n").split("\t")
            if len(column) < 3:
                continue
            protein_id, field, raw_value = column[0], column[1].lower(), column[2].strip()
            if "go" in field or field == "name":
                continue
            if "product" in field:
                result.setdefault(protein_id, []).append(raw_value)
            else:
                value = raw_value.split(":", 1)[1] if ":" in raw_value else raw_value
                result.setdefault(protein_id, []).append(value)
    return result


def parse_go(path: Optional[Path]) -> Dict[str, Set[str]]:
    """Extract GO terms from an iprscan-format annotations file. GO lines
    look like: 'protein_id\\tgo_function\\ttranslation ... |0003743||IEA'.
    """
    result: Dict[str, Set[str]] = {}
    if path is None:
        return result
    with open(path) as fh:
        for line in fh:
            column = line.rstrip("\n").split("\t")
            if len(column) < 3 or "go" not in column[1].lower():
                continue
            go_field = column[2].split("|")
            if len(go_field) < 2 or not go_field[1]:
                continue
            result.setdefault(column[0], set()).add(f"GO:{go_field[1]}")
    return result


def parse_phobius(path: Optional[Path]) -> Tuple[Set[str], Dict[str, int]]:
    """Parse funannotate/phobius native output:
    'SEQUENCE_ID\\tTM\\tSP\\tPrediction'. Returns (secreted_ids,
    {protein_id: n_transmembrane_domains}).
    """
    secreted: Set[str] = set()
    transmembrane: Dict[str, int] = {}
    if path is None:
        return secreted, transmembrane
    with open(path) as fh:
        for line in fh:
            if "PREDICTION" in line.upper():
                continue
            column = line.split()
            if len(column) < 3:
                continue
            protein_id, tm, sp = column[0], column[1], column[2]
            if not tm.isdigit():
                continue
            n_tm = int(tm)
            if sp == "Y" and n_tm == 0:
                secreted.add(protein_id)
            elif n_tm > 0:
                transmembrane[protein_id] = n_tm
    return secreted, transmembrane


def parse_signalp(path: Optional[Path]) -> Set[str]:
    """Parse SignalP 6 native output: comment lines start with '#', data
    lines are 'protein_id\\tPrediction\\t...'. A protein is secreted when
    Prediction starts with 'SP'.
    """
    secreted: Set[str] = set()
    if path is None:
        return secreted
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            column = line.split("\t")
            if len(column) >= 2 and column[1].startswith("SP"):
                secreted.add(column[0])
    return secreted


def parse_effectorp3(path: Optional[Path]) -> Dict[str, float]:
    """Parse EffectorP 3 output:
    '# Identifier\\tCytoplasmic effector\\tApoplastic effector\\tNon-effector\\tPrediction'
    'SS109918_000002-T1 gene=... seq_id=... type=cds\\t-\\t-\\tY (0.984)\\tNon-effector'
    Returns {protein_id: probability} for anything not predicted
    'Non-effector'.
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


def _clean_antismash_id(raw_id: str) -> str:
    """bagRNA's antiSMASH SMCOG lines carry ids like
    'SS109918_SS109918_000005-T1.cds' (locus prefix duplicated, '.cds'
    suffix); the cluster-hit lines are already clean
    ('SS109918_000002-T1'). Normalize both to the plain protein id.
    """
    protein_id = raw_id[:-4] if raw_id.endswith(".cds") else raw_id
    prefix, sep, rest = protein_id.partition("_")
    if sep and rest.startswith(prefix + "_"):
        protein_id = rest
    return protein_id


def parse_antismash(paths: Iterable[Optional[Path]]) -> Tuple[Set[str], Set[str]]:
    """Parse funannotate-style 'annotations.antismash.txt' /
    'annotations.antismash.clusters.txt' files. Returns
    (protein_ids_in_a_cluster, protein_ids_with_an_smcog_hit).
    """
    cluster_hits: Set[str] = set()
    smcog_hits: Set[str] = set()
    for path in paths:
        if path is None:
            continue
        with open(path) as fh:
            for line in fh:
                column = line.rstrip("\n").split("\t")
                if len(column) < 3:
                    continue
                protein_id = _clean_antismash_id(column[0])
                note = column[2]
                if "SMCOG" in note:
                    smcog_hits.add(protein_id)
                elif "cluster" in note.lower():
                    cluster_hits.add(protein_id)
    return cluster_hits, smcog_hits