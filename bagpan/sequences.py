"""Minimal, stdlib-only FASTA/sequence utilities: parsing a genome assembly
and splicing CDS nucleotide sequences out of it using GFF3-derived exon/CDS
coordinates (bagpan.genes). Used by dN/dS (bagpan.dnds), the one bagpan
analysis that needs nucleotide sequences rather than the protein FASTA
everything else works from.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def parse_fasta(path: "str | Path") -> Dict[str, str]:
    """Parses a (genome) FASTA file into {sequence_id: sequence}."""
    sequences: Dict[str, str] = {}
    name = None
    chunks: List[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(chunks)
                name = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            sequences[name] = "".join(chunks)
    return sequences


def extract_cds_sequence(
    genome: Dict[str, str], contig: str, intervals: List[Tuple[int, int]], strand: str
) -> str:
    """Splices and concatenates CDS nucleotide sequence from the genome,
    given 1-based inclusive (start, end) GFF3 intervals (any order) and
    strand ('+'/'-'). Intervals are sorted by genomic position and
    concatenated in transcription order (reversed for '-' strand genes)
    before reverse-complementing the whole spliced sequence.
    """
    if contig not in genome:
        raise KeyError(f"contig {contig!r} not found in genome FASTA")
    sequence = genome[contig]
    ordered = sorted(intervals, key=lambda iv: iv[0])
    pieces = [sequence[start - 1:end] for start, end in ordered]
    spliced = "".join(pieces)
    return reverse_complement(spliced) if strand == "-" else spliced
