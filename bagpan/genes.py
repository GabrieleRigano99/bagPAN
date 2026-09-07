"""Gene-level structure from bagRNA's final annotated GFF3: which transcripts
belong to which gene, and which transcript is the "representative" one
(longest summed CDS) - used to stop isoform/transcript proliferation from
inflating orthogroup sizes and gene counts (a real bug: a single gene with 3
annotated transcripts was being counted 3x).

Falls back to a filename heuristic (strip a trailing '-T\\d+') when a
species has no GFF3 available, so callers never have to special-case a
missing file.
"""

from __future__ import annotations

import dataclasses
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, Tuple

from bagpan.discover import SpeciesFiles

TRANSCRIPT_FEATURES = {"mRNA", "ncRNA", "tRNA", "rRNA"}
_TRANSCRIPT_SUFFIX_RE = re.compile(r"-T\d+$")


@dataclasses.dataclass
class GeneAnnotation:
    gene_of: Dict[str, str]  # transcript/protein id -> gene id
    representative_of_gene: Dict[str, str]  # gene id -> representative transcript id
    gene_position: Dict[str, Tuple[str, int, str]]  # gene id -> (contig, start, strand)
    from_gff3: bool


def _parse_attributes(field: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for part in field.strip().rstrip(";").split(";"):
        if not part or "=" not in part:
            continue
        key, _, value = part.partition("=")
        attrs[key] = value
    return attrs


def _parse_gff3(path: Path):
    gene_span: Dict[str, Tuple[str, int, int, str]] = {}
    transcript_gene: Dict[str, str] = {}
    transcript_cds_length: Dict[str, int] = defaultdict(int)

    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            contig, _source, feature, start, end, _score, strand, _frame, attr_field = cols
            attrs = _parse_attributes(attr_field)

            if feature == "gene":
                gene_id = attrs.get("ID")
                if gene_id:
                    gene_span[gene_id] = (contig, int(start), int(end), strand)
            elif feature in TRANSCRIPT_FEATURES:
                transcript_id, parent = attrs.get("ID"), attrs.get("Parent")
                if transcript_id and parent:
                    transcript_gene[transcript_id] = parent
            elif feature == "CDS":
                parent = attrs.get("Parent")
                if parent:
                    transcript_cds_length[parent] += int(end) - int(start) + 1

    return gene_span, transcript_gene, transcript_cds_length


def _from_gff3(path: Path) -> GeneAnnotation:
    gene_span, transcript_gene, transcript_cds_length = _parse_gff3(path)

    transcripts_by_gene: Dict[str, list] = defaultdict(list)
    for transcript_id, gene_id in transcript_gene.items():
        transcripts_by_gene[gene_id].append(transcript_id)

    representative_of_gene = {
        gene_id: max(sorted(transcripts), key=lambda t: transcript_cds_length.get(t, 0))
        for gene_id, transcripts in transcripts_by_gene.items()
    }
    gene_position = {
        gene_id: (contig, start, strand) for gene_id, (contig, start, _end, strand) in gene_span.items()
    }
    return GeneAnnotation(
        gene_of=dict(transcript_gene),
        representative_of_gene=representative_of_gene,
        gene_position=gene_position,
        from_gff3=True,
    )


def _fasta_headers(path: Path) -> Iterable[str]:
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                yield line[1:].strip().split()[0]


def _from_protein_ids_fallback(protein_ids: Iterable[str]) -> GeneAnnotation:
    gene_of: Dict[str, str] = {}
    transcripts_by_gene: Dict[str, list] = defaultdict(list)
    for protein_id in protein_ids:
        gene_id = _TRANSCRIPT_SUFFIX_RE.sub("", protein_id)
        gene_of[protein_id] = gene_id
        transcripts_by_gene[gene_id].append(protein_id)

    representative_of_gene = {
        gene_id: sorted(transcripts)[0] for gene_id, transcripts in transcripts_by_gene.items()
    }
    return GeneAnnotation(
        gene_of=gene_of,
        representative_of_gene=representative_of_gene,
        gene_position={},
        from_gff3=False,
    )


def build_gene_annotation(files: SpeciesFiles) -> GeneAnnotation:
    """The single entry point callers should use: parses the species' GFF3
    when available, otherwise falls back to a '-T\\d+'-stripping heuristic
    over its proteome headers.
    """
    if files.gff3 is not None:
        return _from_gff3(files.gff3)
    print(
        f"WARNING: no GFF3 found for {files.name!r} - falling back to a "
        "'-Tn' suffix heuristic to group transcripts into genes"
    )
    return _from_protein_ids_fallback(_fasta_headers(files.proteome))