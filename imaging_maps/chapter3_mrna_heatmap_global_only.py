#!/usr/bin/env python3
"""HER2 smFISH Eu-RGB exploratory maps using ONE shared global display scale.

Inputs expected under root:
  231/HER2-eu.tif
  474/eu.tif
  SKBR3/eu.tif

Outputs:
  - B_shared_scale_only_maps.png
  - B_shared_scale_tile_means.png
  - B_three_cell_global_comparison.png
  - analysis_summary.csv
  - METHODS_AND_LIMITATIONS.txt

Important:
This script removes the within-image normalization panel on purpose.
It keeps ONLY a shared-scale display across the three cell lines.
Results are exploratory visualizations from RGB-rendered exports, not
absolute transcript quantification.
"""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import tifffile
from scipy.ndimage import gaussian_filter

SAMPLES = (
    ('MDA-MB-231', '231/HER2-eu.tif'),
    ('BT474', '474/eu.tif'),
    ('SKBR3', 'SKBR3/eu.tif'),
)


def read_rgb(path: Path) -> np.ndarray:
    image = tifffile.imread(path)
    if image.ndim != 3 or image.shape[-1] not in (3, 4):
        raise ValueError(f'{path}: expected RGB/RGBA (H,W,3/4), received {image.shape}')
    if image.dtype != np.uint8:
        raise ValueError(f'{path}: expected 8-bit rendered RGB export, received {image.dtype}')
    return image[..., :3]


def local_signal(red: np.ndarray, small_sigma: float, background_sigma: float) -> np.ndarray:
    r = red.astype(np.float32)
    return np.maximum(
        gaussian_filter(r, small_sigma) - gaussian_filter(r, background_sigma),
        0.0,
    ).astype(np.float32)


def robust_high(x: np.ndarray, pct: float = 99.5) -> float:
    step = max(1, int(np.ceil(np.sqrt(x.size / 300_000))))
    val = float(np.percentile(x[::step, ::step], pct))
    return max(val, 1e-6)


def tile_mean(x: np.ndarray, tile: int) -> np.ndarray:
    h, w = x.shape
    return np.asarray([
        [float(x[y:y+tile, j:j+tile].mean()) for j in range(0, w, tile)]
        for y in range(0, h, tile)
    ], dtype=np.float32)


def plot_shared_maps(rows, filename: Path):
    fig, axes = plt.subplots(3, 2, figsize=(11.5, 13), constrained_layout=True)
    last = None
    for i, s in enumerate(rows):
        axes[i, 0].imshow(s['rgb'])
        axes[i, 0].set_ylabel(s['name'], fontsize=13, fontweight='bold')
        last = axes[i, 1].imshow(s['shared'], cmap='inferno', vmin=0, vmax=1)
        for a in axes[i]:
            a.set_xticks([])
            a.set_yticks([])
    axes[0, 0].set_title('RGB export', fontsize=12)
    axes[0, 1].set_title('Shared-scale HER2-associated signal map', fontsize=12)
    fig.colorbar(last, ax=axes[:, 1], shrink=0.72, pad=0.02,
                 label='Shared normalized display contrast (0–1)')
    fig.suptitle('Chapter 3 | HER2-associated Eu RGB signal mapping (shared scale only)', fontsize=17)
    fig.supxlabel('All three images use one common display scale. Exploratory display only.', fontsize=10)
    fig.savefig(filename, dpi=190, bbox_inches='tight')
    plt.close(fig)


def plot_tile_means(rows, filename: Path, tile_size: int):
    fig, axes = plt.subplots(3, 1, figsize=(6.6, 13), constrained_layout=True)
    last = None
    for i, s in enumerate(rows):
        tiles = tile_mean(s['shared'], tile_size)
        last = axes[i].imshow(tiles, cmap='inferno', vmin=0, vmax=1,
                              interpolation='nearest', aspect='equal')
        axes[i].set_title(f"{s['name']} | shared-scale tile mean", fontsize=12)
        axes[i].set_ylabel('Tile row')
        axes[i].set_xlabel('Tile column')
    fig.colorbar(last, ax=axes, shrink=0.72, pad=0.02,
                 label='Mean shared normalized display contrast')
    fig.suptitle('Chapter 3 | Tile-level summary heatmaps', fontsize=17)
    fig.savefig(filename, dpi=190, bbox_inches='tight')
    plt.close(fig)


def plot_combined(rows, filename: Path, tile_size: int):
    fig, axes = plt.subplots(3, 3, figsize=(14.5, 11.5), constrained_layout=True)
    last = None
    col_titles = ('RGB export', 'B: shared-scale map', 'B: shared-scale tile mean')
    for i, s in enumerate(rows):
        axes[i, 0].imshow(s['rgb'])
        axes[i, 0].set_ylabel(s['name'], fontsize=13, fontweight='bold')
        last = axes[i, 1].imshow(s['shared'], cmap='inferno', vmin=0, vmax=1)
        axes[i, 2].imshow(tile_mean(s['shared'], tile_size), cmap='inferno',
                          vmin=0, vmax=1, interpolation='nearest', aspect='equal')
        axes[i, 0].set_xticks([]); axes[i, 0].set_yticks([])
        axes[i, 1].set_xticks([]); axes[i, 1].set_yticks([])
    for a, title in zip(axes[0], col_titles):
        a.set_title(title, fontsize=12)
    for i in range(3):
        axes[i, 2].set_xlabel('Tile column')
        axes[i, 2].set_ylabel('Tile row')
    fig.colorbar(last, ax=axes[:, 1:], shrink=0.74, pad=0.02,
                 label='Shared normalized display contrast')
    fig.suptitle('Chapter 3 | HER2-associated Eu RGB signal mapping', fontsize=18)
    fig.supxlabel('Within-image normalization removed. One common display scale is used for all three cell lines.', fontsize=10)
    fig.savefig(filename, dpi=190, bbox_inches='tight')
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--tile', type=int, default=128)
    p.add_argument('--small-sigma', type=float, default=2.0)
    p.add_argument('--background-sigma', type=float, default=20.0)
    p.add_argument('--percentile', type=float, default=99.5)
    args = p.parse_args()

    if args.tile < 1 or args.small_sigma <= 0 or args.background_sigma <= args.small_sigma:
        p.error('Require tile>=1 and 0<small-sigma<background-sigma')
    if not 50 < args.percentile < 100:
        p.error('Require 50 < percentile < 100')

    rows = []
    for name, rel in SAMPLES:
        path = args.root / rel
        if not path.is_file():
            p.error(f'Missing image: {path}')
        rgb = read_rgb(path)
        raw = local_signal(rgb[..., 0], args.small_sigma, args.background_sigma)
        rows.append({'name': name, 'filename': rel, 'rgb': rgb, 'raw': raw})

    sampled = []
    for s in rows:
        step = max(1, int(np.ceil(np.sqrt(s['raw'].size / 300_000))))
        sampled.append(s['raw'][::step, ::step].ravel())
    global_scale = max(float(np.percentile(np.concatenate(sampled), args.percentile)), 1e-6)

    for s in rows:
        s['shared'] = np.clip(s['raw'] / global_scale, 0, 1)
        s['tile_mean_mean'] = float(tile_mean(s['shared'], args.tile).mean())
        s['shared_mean'] = float(s['shared'].mean())
        s['shared_p95'] = float(np.percentile(s['shared'], 95.0))
        s['fraction_above_half'] = float((s['raw'] >= global_scale / 2).mean())

    args.out.mkdir(parents=True, exist_ok=True)
    plot_shared_maps(rows, args.out / 'B_shared_scale_only_maps.png')
    plot_tile_means(rows, args.out / 'B_shared_scale_tile_means.png', args.tile)
    plot_combined(rows, args.out / 'B_three_cell_global_comparison.png', args.tile)

    with (args.out / 'analysis_summary.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow((
            'cell_line', 'source_tiff', 'height_px', 'width_px',
            'shared_scale_percentile', 'shared_scale_value',
            'mean_local_contrast_raw', 'mean_shared_normalized',
            'p95_shared_normalized', 'tile_mean_of_shared_map',
            'fraction_above_half_shared_scale'
        ))
        for s in rows:
            writer.writerow((
                s['name'], s['filename'], s['rgb'].shape[0], s['rgb'].shape[1],
                args.percentile, f'{global_scale:.6f}',
                f"{s['raw'].mean():.6f}", f"{s['shared_mean']:.6f}",
                f"{s['shared_p95']:.6f}", f"{s['tile_mean_mean']:.6f}",
                f"{s['fraction_above_half']:.6f}"
            ))

    with (args.out / 'METHODS_AND_LIMITATIONS.txt').open('w', encoding='utf-8') as f:
        f.write(
            'Exploratory analysis of RGB-rendered Eu images.\n'
            'Within-image normalization panels were intentionally removed.\n'
            'Processing: red channel -> local contrast by Gaussian small sigma minus Gaussian background sigma -> clip negatives.\n'
            f'Parameters: sigma_small={args.small_sigma}; sigma_background={args.background_sigma}; percentile={args.percentile}; tile_px={args.tile}.\n'
            f'One shared percentile-based scale was used across all three images; shared_scale={global_scale:.6f}.\n'
            'This is still a display-level comparison of exported RGB images only. It does not correct for exposure, LUT, gain, optics, sample preparation, or field differences.\n'
            'Therefore results should be described as HER2-associated Eu RGB signal maps, not absolute mRNA abundance.\n'
        )

    print('Completed:', args.out)
    print('Files:')
    print('  B_shared_scale_only_maps.png')
    print('  B_shared_scale_tile_means.png')
    print('  B_three_cell_global_comparison.png')
    print('  analysis_summary.csv')
    print('  METHODS_AND_LIMITATIONS.txt')
    print('Shared scale upper bound:', round(global_scale, 3))
    for s in rows:
        print(f"  {s['name']}: mean_shared={s['shared_mean']:.4f}, tile_mean={s['tile_mean_mean']:.4f}, frac>=half={s['fraction_above_half']:.4f}")
    print('WARNING: this shared-scale RGB mapping is exploratory and not absolute transcript quantification.')


if __name__ == '__main__':
    main()
