"""Locate the proteome and annotation files bagPAN needs inside a completed
bagRNA output directory, and auto-detect the isolate's locus-tag prefix from
the proteome's own fasta headers.

bagRNA no longer uses `funannotate annotate`; its functional-annotation
stage is now modules/annotate_functional.nf (bin/merge_functional_annotations.py),
which merges every tool's output into one gene-level
`functional_annotation/functional_annotation.tsv` plus a final
`functional_annotation/annotated.gff3`. The final protein FASTA comes from
AGAT (`structural_annotation/final_proteins.faa`), not funannotate.
"""

from __future__ import annotations

import dataclasses
from pathlib import Path
from typing import Optional


class DiscoveryError(RuntimeError):
    pass


@dataclasses.dataclass
class SpeciesFiles:
    name: str
    root: Path
    proteome: Path
    locus_prefix: str
    functional_annotation_tsv: Path
    gff3: Optional[Path] = None
    antismash_regions_dir: Optional[Path] = None
    effectorp3: Optional[Path] = None


def _detect_locus_prefix(proteome: Path) -> str:
    with open(proteome) as fh:
        for line in fh:
            if line.startswith(">"):
                header = line[1:].strip().split()[0]
                if "_" not in header:
                    raise DiscoveryError(
                        f"Cannot derive a locus-tag prefix from protein header {header!r} "
                        f"in {proteome} (expected LOCUSPREFIX_geneID-Txx format)"
                    )
                return header.split("_", 1)[0]
    raise DiscoveryError(f"No fasta records found in {proteome}")


def discover_species(name: str, bagrna_outdir: "str | Path") -> SpeciesFiles:
    """Locate the proteome and annotation files bagPAN needs for one species,
    inside a completed bagRNA output directory.
    """
    root = Path(bagrna_outdir).resolve()
    if not root.is_dir():
        raise DiscoveryError(f"bagRNA output directory not found: {root}")

    func_root = root / "functional_annotation"
    struct_root = root / "structural_annotation"

    proteome = struct_root / "final_proteins.faa"
    if not proteome.is_file():
        raise DiscoveryError(
            f"No final_proteins.faa found at {proteome} for species {name!r}. "
            "Has structural annotation (AGAT_EXTRACT_PROTEINS) finished for this run?"
        )
    locus_prefix = _detect_locus_prefix(proteome)

    functional_annotation_tsv = func_root / "functional_annotation.tsv"
    if not functional_annotation_tsv.is_file():
        raise DiscoveryError(
            f"No functional_annotation.tsv found at {functional_annotation_tsv} for species "
            f"{name!r}. Has functional annotation (ANNOTATE_FUNCTIONAL) finished for this run?"
        )

    gff3 = func_root / "annotated.gff3"
    if not gff3.is_file():
        gff3 = None

    effectorp3 = func_root / "effectorp3_output.txt"
    if not effectorp3.is_file():
        effectorp3 = None

    antismash_regions_dir = func_root / "antismash_output"
    if not antismash_regions_dir.is_dir() or not any(antismash_regions_dir.glob("*.region*.gbk")):
        antismash_regions_dir = None

    return SpeciesFiles(
        name=name,
        root=root,
        proteome=proteome,
        locus_prefix=locus_prefix,
        functional_annotation_tsv=functional_annotation_tsv,
        gff3=gff3,
        antismash_regions_dir=antismash_regions_dir,
        effectorp3=effectorp3,
    )
