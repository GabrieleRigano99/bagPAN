from bagpan.sequences import extract_cds_sequence, parse_fasta, reverse_complement


def test_reverse_complement():
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("AAGG") == "CCTT"
    assert reverse_complement("ATGN") == "NCAT"


def test_parse_fasta(tmp_path):
    path = tmp_path / "genome.fa"
    path.write_text(">ctg1 some description\nACGTACGT\nACGT\n>ctg2\nTTTTGGGG\n")
    sequences = parse_fasta(path)
    assert sequences == {"ctg1": "ACGTACGTACGT", "ctg2": "TTTTGGGG"}


def test_extract_cds_sequence_plus_strand_single_exon():
    genome = {"ctg1": "AAAATGCCCTAAAAAA"}
    #                   123456789012345 6   (1-based)
    # positions 4-12 = "ATGCCCTAA"
    seq = extract_cds_sequence(genome, "ctg1", [(4, 12)], "+")
    assert seq == "ATGCCCTAA"


def test_extract_cds_sequence_plus_strand_multi_exon_splicing():
    genome = {"ctg1": "ATGCCC" + "GGGGG" + "TAAAAA"}
    # exon1 = positions 1-6 ("ATGCCC"), intron skipped, exon2 = positions 12-17 ("TAAAAA")
    seq = extract_cds_sequence(genome, "ctg1", [(12, 17), (1, 6)], "+")
    assert seq == "ATGCCCTAAAAA"  # intervals sorted by position regardless of input order


def test_extract_cds_sequence_minus_strand_reverse_complements_and_reorders():
    # gene transcribed right-to-left: genomic exon2 (upstream in coords) comes
    # SECOND in transcription order, exon1 (downstream in coords) comes first
    genome = {"ctg1": "TTACAT" + "GGGGG" + "CCCGGG"}
    #                   1-6        7-11     12-17
    # On the minus strand, transcription starts at the high-coordinate end:
    # exon at 12-17 ("CCCGGG") is 5' in transcription order, exon at 1-6 ("TTACAT") is 3'.
    seq = extract_cds_sequence(genome, "ctg1", [(1, 6), (12, 17)], "-")
    # genomic (sorted by position): "TTACAT" + "CCCGGG" -> reverse complement of the whole thing
    from bagpan.sequences import reverse_complement as rc
    assert seq == rc("TTACAT" + "CCCGGG")


def test_extract_cds_sequence_unknown_contig_raises():
    import pytest

    with pytest.raises(KeyError):
        extract_cds_sequence({}, "missing", [(1, 3)], "+")
