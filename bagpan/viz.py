"""Static, dependency-free visualization: hand-rolled SVG (no matplotlib) for
a presence/absence matrix, per-category enrichment bars, and the pangenome
accumulation curve, assembled into a single local report.html. Brings bagpan
to parity with the baseline visual deliverable every mature comparator ships
(Roary's roary_plots.py, PPanGGOLiN's graph, panX's browser) without
reintroducing the heavy plotting deps bagpan otherwise avoids.
"""

from __future__ import annotations

from html import escape as _esc
from pathlib import Path
from typing import Dict, Optional, Tuple

from bagpan.orthogroups import OrthogroupSet

_CATEGORY_COLOR = {"Core": "#d62728", "Accessory": "#ff7f0e", "Singleton": "#f2c744"}
_CATEGORY_ORDER = {"Core": 0, "Accessory": 1, "Singleton": 2}
_PRESENT_COLOR = "#2c7fb8"
_ABSENT_COLOR = "#e8e8e8"


def presence_absence_matrix_svg(
    orthogroup_set: OrthogroupSet, cell_size: int = 4, max_orthogroups: int = 3000
) -> str:
    species = list(orthogroup_set.species_names)
    clusters = sorted(
        orthogroup_set.proteins,
        key=lambda c: (
            _CATEGORY_ORDER[orthogroup_set.category[c]],
            -len(orthogroup_set.species_present[c]),
            c,
        ),
    )
    truncated = len(clusters) > max_orthogroups
    if truncated:
        clusters = clusters[:max_orthogroups]

    margin_left = 10 + max((len(s) for s in species), default=0) * 6
    margin_top = 34
    width = margin_left + len(clusters) * cell_size + 20
    height = margin_top + len(species) * cell_size + 46

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'font-family="sans-serif" font-size="10">'
    ]
    for i, sp in enumerate(species):
        y = margin_top + i * cell_size + cell_size * 0.85
        svg.append(f'<text x="2" y="{y:.1f}">{_esc(sp)}</text>')

    x = margin_left
    for cluster in clusters:
        color = _CATEGORY_COLOR[orthogroup_set.category[cluster]]
        svg.append(f'<rect x="{x}" y="{margin_top - 8}" width="{cell_size}" height="6" fill="{color}"/>')
        present = orthogroup_set.species_present[cluster]
        for i, sp in enumerate(species):
            fill = _PRESENT_COLOR if sp in present else _ABSENT_COLOR
            y = margin_top + i * cell_size
            svg.append(f'<rect x="{x}" y="{y}" width="{cell_size}" height="{cell_size}" fill="{fill}"/>')
        x += cell_size

    legend_y = margin_top + len(species) * cell_size + 16
    lx = margin_left
    for label in ("Core", "Accessory", "Singleton"):
        svg.append(f'<rect x="{lx}" y="{legend_y}" width="10" height="10" fill="{_CATEGORY_COLOR[label]}"/>')
        svg.append(f'<text x="{lx + 14}" y="{legend_y + 9}">{label}</text>')
        lx += 90
    if truncated:
        svg.append(
            f'<text x="{margin_left}" y="{legend_y + 22}">'
            f"(showing first {max_orthogroups} of {len(orthogroup_set.proteins)} orthogroups)</text>"
        )
    svg.append("</svg>")
    return "\n".join(svg)


def category_bar_svg(
    tallies: Dict[str, Dict[str, Tuple[int, int]]], width: int = 560, row_height: int = 16
) -> str:
    categories = list(tallies)
    classes = ("Core", "Accessory", "Singleton")
    block_height = 24 + len(classes) * row_height
    height = 20 + len(categories) * block_height
    max_bar_width = width - 200

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'font-family="sans-serif" font-size="11">'
    ]
    y = 10
    for category in categories:
        svg.append(f'<text x="10" y="{y + 12}" font-weight="bold">{_esc(category)}</text>')
        for j, cls in enumerate(classes):
            hit, total = tallies[category].get(cls, (0, 0))
            pct = (hit / total * 100) if total else 0.0
            bar_width = max_bar_width * (pct / 100)
            by = y + 18 + j * row_height
            svg.append(f'<text x="20" y="{by + 10}">{cls}</text>')
            svg.append(
                f'<rect x="90" y="{by}" width="{bar_width:.1f}" height="12" fill="{_CATEGORY_COLOR[cls]}"/>'
            )
            svg.append(f'<text x="{90 + bar_width + 4:.1f}" y="{by + 10}">{pct:.1f}% ({hit}/{total})</text>')
        y += block_height
    svg.append("</svg>")
    return "\n".join(svg)


def accumulation_curve_svg(
    curve: Dict[int, Dict[str, float]],
    core_fit: Optional[Tuple[float, float]] = None,
    pan_fit: Optional[Tuple[float, float]] = None,
    width: int = 480,
    height: int = 320,
) -> str:
    ks = sorted(curve)
    margin = 42
    plot_w, plot_h = width - 2 * margin, height - 2 * margin
    max_y = max((curve[k]["pan_mean"] for k in ks), default=1) * 1.1 or 1
    max_x = max(ks) if ks else 1

    def sx(k: int) -> float:
        return margin + (k - 1) / (max_x - 1 if max_x > 1 else 1) * plot_w

    def sy(v: float) -> float:
        return height - margin - (v / max_y) * plot_h

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'font-family="sans-serif" font-size="10">'
    ]
    svg.append(f'<line x1="{margin}" y1="{height - margin}" x2="{width - margin}" y2="{height - margin}" stroke="black"/>')
    svg.append(f'<line x1="{margin}" y1="{margin}" x2="{margin}" y2="{height - margin}" stroke="black"/>')
    svg.append(f'<text x="{width / 2 - 60:.1f}" y="{height - 8}">Genomes sampled (k)</text>')

    for k in ks:
        svg.append(f'<circle cx="{sx(k):.1f}" cy="{sy(curve[k]["pan_mean"]):.1f}" r="3" fill="#ff7f0e"/>')
        svg.append(f'<circle cx="{sx(k):.1f}" cy="{sy(curve[k]["core_mean"]):.1f}" r="3" fill="#d62728"/>')

    for fit, color in ((pan_fit, "#ff7f0e"), (core_fit, "#d62728")):
        if fit is None or not ks:
            continue
        kappa, gamma = fit
        points = " ".join(f"{sx(k):.1f},{sy(kappa * (k ** gamma)):.1f}" for k in ks)
        svg.append(f'<polyline points="{points}" fill="none" stroke="{color}" stroke-dasharray="3,2"/>')

    svg.append(f'<circle cx="{width - 130}" cy="16" r="3" fill="#ff7f0e"/><text x="{width - 122}" y="20">Pan-genome</text>')
    svg.append(f'<circle cx="{width - 130}" cy="32" r="3" fill="#d62728"/><text x="{width - 122}" y="36">Core genome</text>')
    svg.append("</svg>")
    return "\n".join(svg)


def _palette(n: int) -> list:
    return [f"hsl({(i * 360 / max(n, 1)) % 360:.0f}, 65%, 55%)" for i in range(n)]


def stacked_bar_svg(
    matrix: Dict[str, Dict[str, int]], bar_width: int = 60, gap: int = 30, plot_height: int = 300
) -> str:
    """Per-species stacked bar chart of gene counts by class (CAZyme family,
    MEROPS class, COG category, BGC type, ...) - funannotate compare's
    CAZy.graph.pdf/COGS.graph.pdf/SM.graph.pdf style, minus matplotlib.
    """
    species = sorted(matrix)
    classes = sorted({c for tally in matrix.values() for c in tally})
    if not species or not classes:
        return '<svg xmlns="http://www.w3.org/2000/svg" width="220" height="40"><text x="6" y="22">(no data)</text></svg>'

    totals = {sp: sum(matrix[sp].values()) for sp in species}
    max_total = max(totals.values()) or 1
    class_color = dict(zip(classes, _palette(len(classes))))

    margin_left = 40
    margin_top = 16
    legend_width = 130
    width = margin_left + len(species) * (bar_width + gap) + legend_width
    height = margin_top + plot_height + 50

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'font-family="sans-serif" font-size="11">'
    ]
    svg.append(
        f'<line x1="{margin_left}" y1="{margin_top}" x2="{margin_left}" '
        f'y2="{margin_top + plot_height}" stroke="black"/>'
    )
    svg.append(
        f'<line x1="{margin_left}" y1="{margin_top + plot_height}" '
        f'x2="{width - legend_width}" y2="{margin_top + plot_height}" stroke="black"/>'
    )

    x = margin_left + gap / 2
    for sp in species:
        y = margin_top + plot_height
        for cls in classes:
            count = matrix[sp].get(cls, 0)
            if count == 0:
                continue
            seg_height = plot_height * (count / max_total)
            y -= seg_height
            svg.append(
                f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_width}" height="{seg_height:.1f}" '
                f'fill="{class_color[cls]}"/>'
            )
        svg.append(
            f'<text x="{x + bar_width / 2:.1f}" y="{margin_top + plot_height + 16}" '
            f'text-anchor="middle">{_esc(sp)}</text>'
        )
        svg.append(
            f'<text x="{x + bar_width / 2:.1f}" y="{margin_top + plot_height + 30}" '
            f'text-anchor="middle">{totals[sp]}</text>'
        )
        x += bar_width + gap

    legend_x = width - legend_width + 10
    legend_y = margin_top
    for cls in classes:
        svg.append(f'<rect x="{legend_x}" y="{legend_y}" width="10" height="10" fill="{class_color[cls]}"/>')
        svg.append(f'<text x="{legend_x + 14}" y="{legend_y + 9}">{_esc(cls)}</text>')
        legend_y += 16
    svg.append("</svg>")
    return "\n".join(svg)


def write_report_html(
    outdir: Path,
    orthogroup_set: OrthogroupSet,
    pangenome_counts: Dict[str, int],
    tallies: Dict[str, Dict[str, Tuple[int, int]]],
    curve: Optional[Dict[int, Dict[str, float]]],
    core_fit: Optional[Tuple[float, float]],
    pan_fit: Optional[Tuple[float, float]],
) -> None:
    matrix_svg = presence_absence_matrix_svg(orthogroup_set)
    bar_svg = category_bar_svg(tallies)
    curve_svg = accumulation_curve_svg(curve, core_fit, pan_fit) if curve else "<p>(not computed)</p>"
    species_list = ", ".join(_esc(s) for s in orthogroup_set.species_names)

    file_links = [
        "pangenome_stats.tsv", "gene_counts.tsv", "orthogroup_classification.tsv",
        "functional_enrichment_summary.tsv", "per_protein_annotations.tsv",
        "pangenome_accumulation.tsv", "pangenome_openness.txt", "run_manifest.json",
    ]
    links_html = "\n".join(f'<a href="{f}">{f}</a>' for f in file_links if (outdir / f).exists())

    html = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>bagPAN report</title>
<style>
body {{ font-family: sans-serif; margin: 24px; color: #222; }}
h2 {{ margin-top: 2em; }}
table {{ border-collapse: collapse; }}
td, th {{ padding: 4px 10px; border: 1px solid #ccc; text-align: left; }}
.files a {{ margin-right: 14px; }}
</style></head>
<body>
<h1>bagPAN report</h1>
<p>{len(orthogroup_set.species_names)} species: {species_list}</p>
<table>
<tr><th>Pangenome category</th><th>Orthogroups</th></tr>
<tr><td>Core</td><td>{pangenome_counts.get("Core", 0)}</td></tr>
<tr><td>Accessory</td><td>{pangenome_counts.get("Accessory", 0)}</td></tr>
<tr><td>Singleton</td><td>{pangenome_counts.get("Singleton", 0)}</td></tr>
</table>
<h2>Presence / absence matrix</h2>
{matrix_svg}
<h2>Functional category enrichment</h2>
{bar_svg}
<h2>Pangenome accumulation curve</h2>
{curve_svg}
<h2>Full output files</h2>
<p class="files">{links_html}</p>
</body></html>
"""
    (outdir / "report.html").write_text(html)