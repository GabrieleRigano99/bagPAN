"""Opt-in BiG-SCAPE integration for secondary-metabolite gene-cluster-family
(GCF) clustering, replacing bagpan's crude presence/absence
"secondary_metabolites" category with actual GCF ids where resolvable.

Mirrors orthofinder_runner.py's Docker-shelling pattern against BiG-SCAPE
2.0's first-party image (ghcr.io/medema-group/big-scape). Unlike the
OrthoFinder integration, this has NOT been validated against a real run in
this codebase (that would need a large Pfam-A.hmm database and real
antiSMASH output) - its SQLite output is parsed defensively: table/column
names are introspected at runtime rather than hardcoded, and a schema this
code doesn't recognize degrades to a clear warning rather than a crash, with
the rest of the bagpan report unaffected.
"""

from __future__ import annotations

import dataclasses
import os
import shutil
import sqlite3
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

from bagpan.discover import SpeciesFiles

DEFAULT_IMAGE = "ghcr.io/medema-group/big-scape:2.0.0-beta.6"


class BigscapeError(RuntimeError):
    pass


@dataclasses.dataclass
class BigscapeResult:
    output_dir: Path
    database: Optional[Path]
    gcf_by_bgc: Dict[str, str]


def stage_bgc_regions(species_files: List[SpeciesFiles], input_dir: Path) -> List[SpeciesFiles]:
    """Copies each species' antiSMASH per-region GenBank files into a single
    flat input dir (species-prefixed to avoid filename collisions - many
    species will independently have e.g. 'region001.gbk'). Returns the
    subset of species_files that actually had region files to contribute.
    """
    input_dir.mkdir(parents=True, exist_ok=True)
    usable = []
    for sf in species_files:
        if sf.antismash_regions_dir is None:
            continue
        region_files = sorted(sf.antismash_regions_dir.glob("*.region*.gbk"))
        if not region_files:
            continue
        for gbk in region_files:
            shutil.copyfile(gbk, input_dir / f"{sf.name}__{gbk.name}")
        usable.append(sf)
    return usable


def _extract_gcf_mapping(db_path: Path) -> Dict[str, str]:
    """Best-effort BGC-id -> gene-cluster-family-id mapping: introspects the
    database's own table/column names rather than assuming a fixed schema,
    since BiG-SCAPE 2.0's exact SQLite layout wasn't verifiable without a
    real run. Returns {} (with a printed warning) if nothing recognizable is
    found - callers should treat that as "no GCF data available", not an
    error.
    """
    mapping: Dict[str, str] = {}
    try:
        conn = sqlite3.connect(str(db_path))
        try:
            cursor = conn.cursor()
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
            tables = [row[0] for row in cursor.fetchall()]
            candidate_tables = [
                t for t in tables if any(k in t.lower() for k in ("family", "gcf", "cluster"))
            ]
            for table in candidate_tables:
                cursor.execute(f"PRAGMA table_info({table})")
                columns = [row[1] for row in cursor.fetchall()]
                id_col = next(
                    (c for c in columns if any(k in c.lower() for k in ("bgc", "region", "gbk"))), None
                )
                family_col = next(
                    (c for c in columns if any(k in c.lower() for k in ("family", "gcf"))), None
                )
                if not (id_col and family_col):
                    continue
                cursor.execute(f"SELECT {id_col}, {family_col} FROM {table}")
                for bgc_id, family_id in cursor.fetchall():
                    if bgc_id is not None and family_id is not None:
                        mapping[str(bgc_id)] = str(family_id)
                if mapping:
                    break
        finally:
            conn.close()
    except sqlite3.Error as exc:
        print(f"WARNING: could not read BiG-SCAPE database {db_path}: {exc}")

    if not mapping:
        print(
            f"WARNING: no recognizable BGC -> gene-cluster-family mapping found in {db_path}. "
            "BiG-SCAPE's schema may differ from what bagpan expects (unvalidated integration) - "
            "secondary-metabolite orthogroups will still get presence/absence, just without a "
            "resolved GCF id."
        )
    return mapping


def run_bigscape(
    species_files: List[SpeciesFiles],
    workdir: Path,
    pfam_hmm_path: "str | Path",
    image: str = DEFAULT_IMAGE,
    threads: int = 4,
) -> BigscapeResult:
    if shutil.which("docker") is None:
        raise BigscapeError("docker is required to run BiG-SCAPE but was not found on $PATH")

    root = workdir / "bigscape"
    input_dir = root / "input"
    output_dir = root / "output"
    output_dir.mkdir(parents=True, exist_ok=True)

    usable = stage_bgc_regions(species_files, input_dir)
    if not usable:
        raise BigscapeError(
            "none of the given species have antiSMASH region .gbk files "
            "(functional_annotation/antismash_output/*.region*.gbk) to run BiG-SCAPE on"
        )

    pfam_hmm_path = Path(pfam_hmm_path).resolve()
    if not pfam_hmm_path.is_file():
        raise BigscapeError(f"--bigscape-pfam path does not exist: {pfam_hmm_path}")
    pfam_dir = pfam_hmm_path.parent

    cmd = [
        "docker", "run", "--rm",
        "-u", f"{os.getuid()}:{os.getgid()}",
        "-v", f"{root}:/home/data",
        "-v", f"{pfam_dir}:/home/pfam:ro",
        image,
        "cluster",
        "-i", "/home/data/input",
        "-o", "/home/data/output",
        "-p", f"/home/pfam/{pfam_hmm_path.name}",
        "-c", str(threads),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise BigscapeError(
            f"BiG-SCAPE Docker run failed (exit {result.returncode}). "
            f"Command: {' '.join(cmd)}\n--- stderr (tail) ---\n{result.stderr[-4000:]}"
        )

    db_hits = sorted(
        list(output_dir.rglob("*.db")) + list(output_dir.rglob("*.sqlite")) + list(output_dir.rglob("*.sqlite3"))
    )
    if not db_hits:
        print(f"WARNING: BiG-SCAPE ran but no SQLite database was found under {output_dir}")
        return BigscapeResult(output_dir=output_dir, database=None, gcf_by_bgc={})

    database = db_hits[-1]
    return BigscapeResult(output_dir=output_dir, database=database, gcf_by_bgc=_extract_gcf_mapping(database))