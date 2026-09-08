"""Pairwise dN/dS (Nei & Gojobori 1986 method) between single-copy orthologs,
computed entirely in pure Python - bagpan's replacement for funannotate
compare's `--run_dnds estimate` mode (which shells out to mafft/trimal/PAML),
matching bagpan's stdlib-only design.

Pipeline: translate both CDS -> align the two proteins (Needleman-Wunsch,
linear gap penalty - a deliberate simplicity/reliability tradeoff over a
substitution-matrix + affine-gap aligner, reasonable here since this only
ever runs on already-orthologous, closely related sequences, not remote
homology search) -> back-translate the alignment to paired codons -> classic
Nei-Gojobori synonymous/nonsynonymous site and difference counting, with a
Jukes-Cantor correction for multiple hits.
"""

from __future__ import annotations

import dataclasses
import itertools
import math
from typing import Dict, List, Optional, Tuple

_BASES = "TCAG"
_CODONS = [a + b + c for a in _BASES for b in _BASES for c in _BASES]
_AMINO_ACIDS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
GENETIC_CODE: Dict[str, str] = dict(zip(_CODONS, _AMINO_ACIDS))
STOP = "*"
_VALID_BASES = set("ACGT")


def translate_codon(codon: str) -> str:
    return GENETIC_CODE.get(codon.upper(), "X")


def align_proteins_nw(seq_a: str, seq_b: str, match: int = 2, mismatch: int = -1, gap: int = -2) -> Tuple[str, str]:
    """Global pairwise alignment (Needleman-Wunsch, linear gap penalty)."""
    n, m = len(seq_a), len(seq_b)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        score[i][0] = score[i - 1][0] + gap
    for j in range(1, m + 1):
        score[0][j] = score[0][j - 1] + gap
    for i in range(1, n + 1):
        row, prev_row = score[i], score[i - 1]
        for j in range(1, m + 1):
            diag = prev_row[j - 1] + (match if seq_a[i - 1] == seq_b[j - 1] else mismatch)
            row[j] = max(diag, prev_row[j] + gap, row[j - 1] + gap)

    aligned_a: List[str] = []
    aligned_b: List[str] = []
    i, j = n, m
    while i > 0 or j > 0:
        diag = score[i - 1][j - 1] + (match if seq_a[i - 1] == seq_b[j - 1] else mismatch) if i > 0 and j > 0 else None
        if i > 0 and j > 0 and score[i][j] == diag:
            aligned_a.append(seq_a[i - 1])
            aligned_b.append(seq_b[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and score[i][j] == score[i - 1][j] + gap:
            aligned_a.append(seq_a[i - 1])
            aligned_b.append("-")
            i -= 1
        else:
            aligned_a.append("-")
            aligned_b.append(seq_b[j - 1])
            j -= 1
    return "".join(reversed(aligned_a)), "".join(reversed(aligned_b))


def codon_align(
    aligned_a: str, codons_a: List[str], aligned_b: str, codons_b: List[str]
) -> Tuple[List[str], List[str]]:
    """Walks a protein-level alignment (with '-' gaps) and each sequence's
    own ordered codon list, returning parallel codon lists for columns where
    BOTH sequences have a real (non-gap, unambiguous) residue - the only
    columns usable for Nei-Gojobori counting.
    """
    out_a: List[str] = []
    out_b: List[str] = []
    ia = ib = 0
    for ra, rb in zip(aligned_a, aligned_b):
        codon_a = codons_a[ia] if ra != "-" else None
        codon_b = codons_b[ib] if rb != "-" else None
        if ra != "-":
            ia += 1
        if rb != "-":
            ib += 1
        if codon_a is None or codon_b is None:
            continue
        if not (set(codon_a) <= _VALID_BASES and set(codon_b) <= _VALID_BASES):
            continue
        out_a.append(codon_a)
        out_b.append(codon_b)
    return out_a, out_b


def _synonymous_fraction(codon: str, position: int) -> float:
    """Fraction of the 3 possible single-nt substitutions at `position` that
    are synonymous (same amino acid). A substitution to a stop codon counts
    as nonsynonymous (the common NG86 convention).
    """
    original_aa = translate_codon(codon)
    syn = 0
    for base in "ACGT":
        if base == codon[position]:
            continue
        mutant_aa = translate_codon(codon[:position] + base + codon[position + 1:])
        if mutant_aa == original_aa:
            syn += 1
    return syn / 3.0


def count_sites(codon: str) -> Tuple[float, float]:
    """(synonymous_sites, nonsynonymous_sites) for one codon, NG86-style."""
    syn = sum(_synonymous_fraction(codon, p) for p in range(3))
    return syn, 3.0 - syn


def count_differences_ng86(codon_a: str, codon_b: str) -> Tuple[float, float]:
    """(synonymous_diffs, nonsynonymous_diffs) between two codons, averaged
    over all shortest mutational pathways connecting them (Nei-Gojobori);
    pathways that pass through a stop codon are excluded.
    """
    diff_positions = [p for p in range(3) if codon_a[p] != codon_b[p]]
    if not diff_positions:
        return 0.0, 0.0

    syn_total = nonsyn_total = 0.0
    n_valid_pathways = 0
    for order in itertools.permutations(diff_positions):
        current = codon_a
        path_syn = path_nonsyn = 0
        valid = True
        for pos in order:
            next_codon = current[:pos] + codon_b[pos] + current[pos + 1:]
            aa_before, aa_after = translate_codon(current), translate_codon(next_codon)
            if aa_after == STOP:
                valid = False
                break
            path_syn += aa_after == aa_before
            path_nonsyn += aa_after != aa_before
            current = next_codon
        if not valid:
            continue
        syn_total += path_syn
        nonsyn_total += path_nonsyn
        n_valid_pathways += 1

    if n_valid_pathways == 0:
        return 0.0, float(len(diff_positions))
    return syn_total / n_valid_pathways, nonsyn_total / n_valid_pathways


def _jukes_cantor_correct(p: float) -> Optional[float]:
    if p >= 0.75:
        return None
    return -0.75 * math.log(1 - (4.0 / 3.0) * p)


@dataclasses.dataclass
class DnDsResult:
    n_codons_compared: int
    syn_sites: float
    nonsyn_sites: float
    syn_diffs: float
    nonsyn_diffs: float
    dS: Optional[float]
    dN: Optional[float]
    omega: Optional[float]


def _codons_and_protein(cds: str) -> Tuple[List[str], str]:
    cds = cds.upper()
    n_codons = len(cds) // 3
    codons = [cds[i:i + 3] for i in range(0, n_codons * 3, 3)]
    protein = "".join(translate_codon(c) for c in codons)
    if protein.endswith(STOP):
        protein = protein[:-1]
        codons = codons[:-1]
    return codons, protein


def pairwise_dn_ds(cds_a: str, cds_b: str) -> DnDsResult:
    """Full pipeline: translate -> align -> codon-align -> NG86 count ->
    Jukes-Cantor-corrected dS/dN/omega for one pair of CDS sequences.
    """
    codons_a, protein_a = _codons_and_protein(cds_a)
    codons_b, protein_b = _codons_and_protein(cds_b)
    aligned_a, aligned_b = align_proteins_nw(protein_a, protein_b)
    paired_a, paired_b = codon_align(aligned_a, codons_a, aligned_b, codons_b)

    syn_sites = nonsyn_sites = 0.0
    syn_diffs = nonsyn_diffs = 0.0
    for ca, cb in zip(paired_a, paired_b):
        sa, na = count_sites(ca)
        sb, nb = count_sites(cb)
        syn_sites += (sa + sb) / 2
        nonsyn_sites += (na + nb) / 2
        sd, nd = count_differences_ng86(ca, cb)
        syn_diffs += sd
        nonsyn_diffs += nd

    pS = (syn_diffs / syn_sites) if syn_sites else None
    pN = (nonsyn_diffs / nonsyn_sites) if nonsyn_sites else None
    dS = _jukes_cantor_correct(pS) if pS is not None else None
    dN = _jukes_cantor_correct(pN) if pN is not None else None
    omega = (dN / dS) if (dN is not None and dS is not None and dS > 0) else None

    return DnDsResult(
        n_codons_compared=len(paired_a),
        syn_sites=syn_sites, nonsyn_sites=nonsyn_sites,
        syn_diffs=syn_diffs, nonsyn_diffs=nonsyn_diffs,
        dS=dS, dN=dN, omega=omega,
    )
