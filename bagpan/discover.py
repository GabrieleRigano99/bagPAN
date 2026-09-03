"""Locate the proteome and annotation files bagPAN needs inside a completed
bagRNA output directory, and auto-detect the isolate's locus-tag prefix from
the proteome's own fasta headers (bagRNA's functional-annotation stage wraps
funannotate, so these paths follow funannotate's output layout).
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
    gff3: Optional[Path] = None
    antismash_regions_dir: Optional[Path] = None
    pfam: Optional[Path] = None
    iprscan: Optional[Path] = None
    dbcan: Optional[Path] = None
    merops: Optional[Path] = None
    genes_products: Optional[Path] = None
    antismash: Optional[Path] = None
    antismash_clusters: Optional[Path] = None
    phobius: Optional[Path] = None
    signalp: Optional[Path] = None
    effectorp3: Optional[Path] = None


def _first_match(patterns, base: Path) -> Optional[Path]:
    if not base.is_dir():
        return None
    for pattern in patterns:
        hits = sorted(base.glob(pattern))
        if hits:
            return hits[0]
    return None


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
    annotate_results = func_root / "functional_annotation" / "annotate_results"
    annotate_misc = func_root / "functional_annotation" / "annotate_misc"

    proteome = _first_match(["*.proteins.fa", "*.proteins.faa"], annotate_results)
    if proteome is None:
        raise DiscoveryError(
            f"No *.proteins.fa found under {annotate_results} for species {name!r}. "
            "Has functional annotation finished for this run?"
        )
    locus_prefix = _detect_locus_prefix(proteome)

    antismash_regions_dir = func_root / "antismash_output"
    if not antismash_regions_dir.is_dir() or not any(antismash_regions_dir.glob("*.region*.gbk")):
        antismash_regions_dir = None

    return SpeciesFiles(
        name=name,
        root=root,
        proteome=proteome,
        locus_prefix=locus_prefix,
        gff3=_first_match(["*.gff3"], annotate_results),
        antismash_regions_dir=antismash_regions_dir,
        pfam=_first_match(["annotations.pfam.txt"], annotate_misc),
        iprscan=_first_match(["annotations.iprscan.txt"], annotate_misc),
        dbcan=_first_match(["annotations.dbCAN.txt", "annotations.dbcan.txt"], annotate_misc),
        merops=_first_match(["annotations.merops.txt"], annotate_misc),
        genes_products=_first_match(["annotations.genes-products.txt"], annotate_misc),
        antismash=_first_match(["annotations.antismash.txt"], annotate_misc),
        antismash_clusters=_first_match(["annotations.antismash.clusters.txt"], annotate_misc),
        phobius=_first_match(["phobius_results.txt"], func_root),
        signalp=_first_match(["signalp_results.txt"], func_root),
        effectorp3=_first_match(["effectorp3_output.txt"], func_root),
    )