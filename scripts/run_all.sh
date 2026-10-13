#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: bash scripts/run_all.sh '/path/to/fig3 0710' '/path/to/Her2 mRNA mapping'" >&2
  exit 2
fi
GRADIENT_ROOT="$1"
MAPPING_ROOT="$2"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
[[ -d "$GRADIENT_ROOT/25%" && -d "$GRADIENT_ROOT/100%" ]] || { echo "Missing gradient directories: $GRADIENT_ROOT" >&2; exit 1; }
[[ -f "$MAPPING_ROOT/231/HER2-eu.tif" && -f "$MAPPING_ROOT/474/eu.tif" && -f "$MAPPING_ROOT/SKBR3/eu.tif" ]] || { echo "Missing one or more expected Eu TIFF: $MAPPING_ROOT" >&2; exit 1; }

mkdir -p "$REPO/outputs" "$REPO/figures"
python "$REPO/imaging_maps/chapter3_mrna_heatmap_global_only.py" \
  --root "$MAPPING_ROOT" --out "$REPO/outputs/mapping" --tile 128
python "$REPO/spot_counts/chapter3_spot_scatter_stats.py" \
  --root "$GRADIENT_ROOT" --out "$REPO/outputs/spot_counts"
python "$REPO/spot_counts/method_robustness.py" \
  --root "$GRADIENT_ROOT" --out "$REPO/outputs/trackmate_robustness"

# Copy reproducible, derived PNG figures into the curated GitHub figures directory.
cp "$REPO/outputs/mapping/B_three_cell_global_comparison.png" \
  "$REPO/figures/her2_shared_scale_map.png"
cp "$REPO/outputs/spot_counts/spots_per_cell_scatter_median.png" \
  "$REPO/figures/trackmate_spots_per_cell.png"
cp "$REPO/outputs/trackmate_robustness/trackmate_quality_sensitivity.png" \
  "$REPO/figures/trackmate_quality_sensitivity.png"

echo
printf '%s\n' '========== Per-cell group statistics =========='
cat "$REPO/outputs/spot_counts/gradient_summary.csv"
echo
printf '%s\n' '========== Exclusion audit =========='
python - "$REPO/outputs/spot_counts/file_audit.csv" <<'PYTHON'
import csv, sys
from collections import Counter
with open(sys.argv[1], encoding='utf-8-sig', newline='') as f:
    rows=list(csv.DictReader(f))
print(dict(Counter(r['status'] for r in rows)))
for r in rows:
    if r['status']!='included':
        print(r['file'], r['reason'])
PYTHON

echo
printf 'Finished. Derived figures: %s/figures/\n' "$REPO"
echo 'Review images and confirm experimental provenance before publishing.'
