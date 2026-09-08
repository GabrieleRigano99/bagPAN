"""bagpan.dnds tests. The core NG86 pipeline (align_proteins_nw + codon_align
+ count_sites + count_differences_ng86 + Jukes-Cantor correction) was
cross-validated during development against biopython's independent
Bio.codonalign NG86 implementation on 15 randomized ~80-codon sequence pairs
(0/15 mismatches beyond floating-point noise) - biopython is not a bagpan
dependency, so that cross-check isn't repeated here; the case below is one
of those validated pairs, captured as a fixed regression value.
"""

import math

from bagpan.dnds import (
    align_proteins_nw,
    codon_align,
    count_differences_ng86,
    count_sites,
    pairwise_dn_ds,
    translate_codon,
)


def test_translate_codon_standard_genetic_code():
    assert translate_codon("ATG") == "M"
    assert translate_codon("TTT") == "F"
    assert translate_codon("TAA") == "*"
    assert translate_codon("TAG") == "*"
    assert translate_codon("TGA") == "*"
    assert translate_codon("GGG") == "G"
    assert translate_codon("atg") == "M"  # case-insensitive
    assert translate_codon("NNN") == "X"  # unknown/ambiguous


def test_align_proteins_nw_identical_sequences():
    aligned_a, aligned_b = align_proteins_nw("MSEQ", "MSEQ")
    assert aligned_a == "MSEQ"
    assert aligned_b == "MSEQ"


def test_align_proteins_nw_single_insertion():
    aligned_a, aligned_b = align_proteins_nw("MSEQAB", "MSEQB")
    assert aligned_a.replace("-", "") == "MSEQAB"
    assert aligned_b.replace("-", "") == "MSEQB"
    assert len(aligned_a) == len(aligned_b)
    # the extra 'A' in seq_a should align against a gap
    assert "-" in aligned_b


def test_codon_align_skips_gapped_columns():
    aligned_a, aligned_b = "MA-Q", "MAEQ"
    codons_a = ["ATG", "GCT", "CAG"]  # M, A, Q  (no codon for the gap)
    codons_b = ["ATG", "GCC", "GAG", "CAA"]  # M, A, E, Q
    paired_a, paired_b = codon_align(aligned_a, codons_a, aligned_b, codons_b)
    assert paired_a == ["ATG", "GCT", "CAG"]
    assert paired_b == ["ATG", "GCC", "CAA"]  # E (gap column in seq_a) skipped


def test_count_sites_known_codon():
    # AAA (Lys): position 2 (3rd base) has exactly 1/3 synonymous substitutions
    # (AAA->AAG is Lys->Lys), positions 0 and 1 are fully nonsynonymous.
    syn, nonsyn = count_sites("AAA")
    assert abs(syn - 1 / 3) < 1e-9
    assert abs(nonsyn - 8 / 3) < 1e-9


def test_count_differences_ng86_single_synonymous_substitution():
    # AAA -> AAG: 3rd-position change, Lys -> Lys
    syn, nonsyn = count_differences_ng86("AAA", "AAG")
    assert syn == 1.0 and nonsyn == 0.0


def test_count_differences_ng86_single_nonsynonymous_substitution():
    # AAA -> ACA: 2nd-position change, Lys -> Thr
    syn, nonsyn = count_differences_ng86("AAA", "ACA")
    assert syn == 0.0 and nonsyn == 1.0


def test_count_differences_ng86_identical_codons():
    assert count_differences_ng86("ATG", "ATG") == (0.0, 0.0)


def test_pairwise_dn_ds_identical_sequences_gives_zero():
    cds = "ATG" + "GCT" * 20 + "TAA"
    result = pairwise_dn_ds(cds, cds)
    assert result.dS == 0.0
    assert result.dN == 0.0
    assert result.omega is None  # 0/0 is undefined, not infinite or zero


def test_pairwise_dn_ds_cross_validated_against_biopython():
    # Captured from a bagpan-vs-biopython(Bio.codonalign, method='NG86')
    # cross-check during development; matched to 4 decimal places.
    base = (
        "TATACGGACGCAGTGTCTCTGTCGATAGTGCGGATCAGTCAGGCCCTTTCAATATTCGGTGATCAGCGA"
        "AACGTGGCTTTTAGGCGGCCTGTCGCACTCAATGGGTCAGGTCCGTTCTTCTTCAGTACCTTTGGGGAG"
        "CAGAGACTTCGAGTCTTCACTCTCGTGCGGGGGATAACACTCCACCTCAGACTCGTGATTGGGCCCGG"
        "ATTCCGCGATGGCACAGGAAGTTCATGCAAGGTCTAA"
    )
    mutant = (
        "TATACGGACGCAGTGACTCTGTCCATAGTGCGGATCAGTCAAGCCCTTTCTATATTCGGTGATCAGCGA"
        "AACGTGGCTTTTAGGCGGCCTGTCGCATTCAATGGGTCAGGTCCGTTCTTCTTCAGTACCTTTGGGGAG"
        "CAGAAACTTCGAGTCTTCACTCTCGTGCGGGGGATAACACTCCACCTCAGACTCGTAATTGGGCCTGG"
        "ATTCCGCGATGGCACAGGAAGTTCAAGCAAGCTCTAA"
    )
    result = pairwise_dn_ds(base, mutant)
    assert result.n_codons_compared == 80
    assert math.isclose(result.dS, 0.0841216928217367, rel_tol=1e-6)
    assert math.isclose(result.dN, 0.028766709292741117, rel_tol=1e-6)
    assert math.isclose(result.omega, 0.34196541139157766, rel_tol=1e-6)


def test_pairwise_dn_ds_saturated_case_returns_none_not_error():
    # 1 codon compared with a synonymous change: pS = 1/(1/3) = 3.0 >= 0.75,
    # Jukes-Cantor correction is undefined there - must return None, not raise.
    result = pairwise_dn_ds("AAA" + "TAA", "AAG" + "TAA")
    assert result.dS is None
    assert result.omega is None
