"""Stage per-species proteomes and run OrthoFinder via Docker."""

from __future__ import annotations

import dataclasses
import os
import shutil
import subprocess
from pathlib import Path
from typing import List, Optional

from bagpan.discover import SpeciesFiles

DEFAULT_IMAGE = "davidemms/orthofinder:2.5.5"


class OrthoFinderError(RuntimeError):
    pass


@dataclasses.dataclass
class OrthoFinderResult:
    orthogroups: Path
    species_tree: Optional[Path]
    statistics_overall: Optional[Path]


def stage_proteomes(species_files: List[SpeciesFiles], staging_dir: Path) -> None:
    staging_dir.mkdir(parents=True, exist_ok=True)
    for sf in species_files:
        shutil.copyfile(sf.proteome, staging_dir / f"{sf.name}.fa")


def run_orthofinder(
    species_files: List[SpeciesFiles],
    workdir: Path,
    image: str = DEFAULT_IMAGE,
    threads: int = 4,
) -> OrthoFinderResult:
    """Stage proteomes under workdir/orthofinder_input and run OrthoFinder in
    a Docker container. Returns the paths to the resulting Orthogroups.txt
    plus, when OrthoFinder produced them, its own rooted species tree and
    overall comparative-genomics statistics - both free byproducts of the
    same run that bagpan surfaces instead of discarding.
    """
    if shutil.which("docker") is None:
        raise OrthoFinderError(
            "docker is required to run OrthoFinder but was not found on $PATH "
            "(use --skip-orthofinder --orthogroups PATH if you already have one)"
        )

    staging_dir = workdir / "orthofinder_input"
    stage_proteomes(species_files, staging_dir)

    cmd = [
        "docker",
        "run",
        "--rm",
        "-u",
        f"{os.getuid()}:{os.getgid()}",
        "-v",
        f"{staging_dir}:/data",
        image,
        "orthofinder",
        "-f",
        "/data",
        "-t",
        str(threads),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise OrthoFinderError(
            f"OrthoFinder Docker run failed (exit {result.returncode}). "
            f"Command: {' '.join(cmd)}\n--- stderr (tail) ---\n{result.stderr[-4000:]}"
        )

    orthogroup_hits = sorted(staging_dir.glob("OrthoFinder/Results_*/Orthogroups/Orthogroups.txt"))
    if not orthogroup_hits:
        raise OrthoFinderError(
            "OrthoFinder ran but no Orthogroups.txt was found under "
            f"{staging_dir}/OrthoFinder/Results_*/Orthogroups/"
        )
    results_dir = orthogroup_hits[-1].parent.parent

    tree_hits = sorted(results_dir.glob("Species_Tree/SpeciesTree_rooted.txt"))
    stats_hits = sorted(results_dir.glob("Comparative_Genomics_Statistics/Statistics_Overall.tsv"))

    return OrthoFinderResult(
        orthogroups=orthogroup_hits[-1],
        species_tree=tree_hits[-1] if tree_hits else None,
        statistics_overall=stats_hits[-1] if stats_hits else None,
    )