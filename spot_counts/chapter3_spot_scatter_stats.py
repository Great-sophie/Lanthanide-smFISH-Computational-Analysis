#!/usr/bin/env python3
"""Chapter 3: per-cell TrackMate spot count scatterplot and descriptive statistics.

Each '<gradient>-<cell_number>.csv' is one cell (experimenter-confirmed).
Only filenames matching their directory's dUTP-N3 percentage are included.
TrackMate metadata occupies the first four rows; blank/footer rows are skipped.
No TIFFs are required. Counts are candidate spots, not verified transcripts.
"""
from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

GROUPS = ("25%", "50%", "75%", "100%")
STEM = re.compile(r"^(25%|50%|75%|100%)-(\d+)$")


def read_spots(path: Path) -> tuple[int, int, int]:
    """Return (n_valid_spots, n_blank_or_footer_rows, n_invalid_rows)."""
    with path.open(encoding="utf-8-sig", errors="replace", newline="") as fh:
        rows = list(csv.reader(fh))
    if len(rows) < 4:
        raise ValueError("Missing four-row TrackMate header/metadata")
    header = [col.strip().upper() for col in rows[0]]
    required = ("ID", "QUALITY", "POSITION_X", "POSITION_Y")
    if not all(col in header for col in required):
        raise ValueError("Not a TrackMate Spot CSV: missing ID/QUALITY/POSITION_X/POSITION_Y")
    indices = {name: header.index(name) for name in required}
    spot_ids: set[str] = set()
    footer = invalid = 0

    for row in rows[4:]:
        def field(name: str) -> str:
            idx = indices[name]
            return row[idx].strip() if idx < len(row) else ""

        sid, quality, x, y = [field(name) for name in required]
        if not any((sid, quality, x, y)):
            # Excel/TrackMate CSVs may have blank rows or footer data in other columns.
            footer += 1
            continue
        try:
            values = (float(quality), float(x), float(y))
            if not sid or not all(math.isfinite(v) for v in values):
                raise ValueError("Empty ID or invalid numeric value")
        except ValueError:
            invalid += 1
            continue
        if sid in spot_ids:
            invalid += 1
            continue
        spot_ids.add(sid)
    return len(spot_ids), footer, invalid


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8-sig") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Single-cell TrackMate spots: scatter, median, sample variance")
    parser.add_argument("--root", required=True, type=Path, help="fig3 0710 folder containing 25%%, 50%%, 75%%, 100%%")
    parser.add_argument("--out", required=True, type=Path, help="Output folder")
    args = parser.parse_args()
    root = args.root.expanduser().resolve()
    out = args.out.expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)

    cells: list[dict] = []
    audit: list[dict] = []
    for group in GROUPS:
        folder = root / group
        if not folder.is_dir():
            parser.error(f"Missing folder: {folder}")
        for file in sorted(folder.glob("*.csv")):
            audit_row = {"folder_group": group, "file": str(file.relative_to(root)), "status": "", "reason": "", "blank_footer_rows": "", "invalid_rows": ""}
            stem_match = STEM.fullmatch(file.stem)
            if stem_match is None:
                audit_row.update(status="excluded", reason="filename does not identify a numbered cell")
            elif stem_match.group(1) != group:
                audit_row.update(status="excluded", reason="filename/folder gradient mismatch")
            else:
                try:
                    count, blanks, invalid = read_spots(file)
                    audit_row["blank_footer_rows"] = blanks
                    audit_row["invalid_rows"] = invalid
                    if invalid:
                        audit_row.update(status="excluded", reason="nonblank malformed or duplicate spot record; review manually")
                    elif count == 0:
                        audit_row.update(status="excluded", reason="zero valid spots; manual review needed")
                    else:
                        audit_row.update(status="included", reason="")
                        cells.append({"gradient": group, "cell_id": file.stem, "spots_per_cell": count, "source_csv": str(file.relative_to(root))})
                except Exception as exc:
                    audit_row.update(status="excluded", reason=f"parse error: {exc}")
            audit.append(audit_row)

    if not cells:
        write_csv(out / "file_audit.csv", list(audit[0]) if audit else ["file", "status", "reason"], audit)
        raise SystemExit("No valid cells; inspect file_audit.csv")

    summary: list[dict] = []
    fig, ax = plt.subplots(figsize=(9, 6))
    for i, group in enumerate(GROUPS, start=1):
        group_cells = [r for r in cells if r["gradient"] == group]
        values = np.array([r["spots_per_cell"] for r in group_cells], dtype=float)
        if not len(values):
            summary.append({"gradient": group, "n_cells": 0, "median": "", "mean": "", "sample_variance": "", "sample_sd": "", "q1": "", "q3": ""})
            continue
        # Deterministic jitter, so the same input produces the same plot.
        offsets = np.linspace(-0.12, 0.12, len(values)) if len(values) > 1 else np.array([0.])
        ax.scatter(i + offsets, values, s=48, alpha=.82, zorder=3)
        median = float(np.median(values))
        ax.hlines(median, i - .23, i + .23, color="black", lw=2.8, zorder=4)
        ax.text(i, float(np.max(values)) + 9, f"n={len(values)}", ha="center", fontsize=10)
        summary.append({
            "gradient": group, "n_cells": len(values), "median": median,
            "mean": float(np.mean(values)),
            "sample_variance": float(np.var(values, ddof=1)) if len(values) > 1 else "",
            "sample_sd": float(np.std(values, ddof=1)) if len(values) > 1 else "",
            "q1": float(np.percentile(values, 25)), "q3": float(np.percentile(values, 75)),
        })
    ax.set_xticks([1, 2, 3, 4], GROUPS)
    ax.set_xlim(.5, 4.5)
    ax.set_xlabel("dUTP-N3 labeling fraction")
    ax.set_ylabel("TrackMate candidate spots per cell")
    ax.set_title("SKBR3 HER2 smFISH — candidate spot counts per cell")
    ax.grid(axis="y", alpha=.25)
    ax.set_axisbelow(True)
    fig.text(.5, .012, "Each dot = one cell (CSV); black line = group median. Cells are not independent biological replicates.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .04, 1, 1))
    fig.savefig(out / "spots_per_cell_scatter_median.png", dpi=240)
    plt.close(fig)

    write_csv(out / "single_cell_counts.csv", ["gradient", "cell_id", "spots_per_cell", "source_csv"], cells)
    write_csv(out / "gradient_summary.csv", ["gradient", "n_cells", "median", "mean", "sample_variance", "sample_sd", "q1", "q3"], summary)
    write_csv(out / "file_audit.csv", ["folder_group", "file", "status", "reason", "blank_footer_rows", "invalid_rows"], audit)

    print(f"Included: {len(cells)} cells; excluded: {len(audit)-len(cells)} CSVs")
    for item in summary:
        if item["n_cells"]:
            print(f"{item['gradient']:4s} n={item['n_cells']:2d}, median={item['median']:.1f}, mean={item['mean']:.2f}, sample variance={item['sample_variance']:.2f}")
    print(f"Saved: {out}")
    print("Note: TrackMate counts are candidate spots, not validated transcript counts; biological batches unknown.")


if __name__ == "__main__":
    main()
