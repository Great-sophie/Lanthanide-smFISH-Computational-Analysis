#!/usr/bin/env python3
"""Audit and visualize three-cell-line HER2 smFISH microscopy exports.

This is a display/channel QC utility, NOT an expression quantifier.
"""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import tifffile

CELL_FILES = [
    ('MDA-MB-231 (231)', '231', ('cy3.tif', 'HER2-eu.tif', 'DAPI+eu.tif')),
    ('BT474 (474)', '474', ('cy3.tif', 'eu.tif', 'dapi+eu.tif')),
    ('SKBR3', 'SKBR3', ('cy3.tif', 'eu.tif', 'dapi+eu.tif')),
]

def read_tiff(path):
    with tifffile.TiffFile(path) as tif:
        a = tif.asarray()
        info = tif.pages[0]
        photometric = str(info.photometric)
    if a.ndim != 3 or a.shape[-1] not in (3, 4):
        raise ValueError(f'Expected RGB/RGBA HxWxC, got {a.shape}')
    if a.shape[-1] == 4:
        a = a[..., :3]
    return a, photometric

def view_channel(rgb, idx):
    a = rgb[..., idx].astype(np.float32)
    low, high = np.percentile(a, [1, 99.5])
    if high <= low:
        high = low + 1.0
    return np.clip((a - low) / (high-low), 0, 1)

def view_rgb(rgb):
    arr = rgb.astype(np.float32)
    if np.issubdtype(rgb.dtype, np.integer):
        return np.clip(arr / np.iinfo(rgb.dtype).max, 0, 1)
    lo, hi = np.percentile(arr, [0, 99.9])
    return np.clip((arr - lo) / max(hi - lo, 1), 0, 1)

def plot_overview(cache, out):
    fig, axes = plt.subplots(3, 3, figsize=(13, 12), constrained_layout=True)
    for row, (label, group, names) in enumerate(CELL_FILES):
        for col, name in enumerate(names):
            ax = axes[row, col]
            key = (group, name)
            if key in cache:
                a, _ = cache[key]
                ax.imshow(view_rgb(a), interpolation='nearest')
                ax.set_title(f'{label}\n{name}\n{a.shape[1]}×{a.shape[0]}', fontsize=11)
            else:
                ax.text(.5, .5, 'MISSING / UNREADABLE', ha='center', va='center', transform=ax.transAxes)
            ax.axis('off')
    fig.suptitle('Chapter 3 | HER2 smFISH imaging inventory (RGB exports)', fontsize=16)
    fig.savefig(out / 'three_cell_line_overview.png', dpi=160)
    plt.close(fig)

def plot_channels(cache, out):
    # One image per RGB channel plus display original, per source image.
    for (group, name), (a, _) in cache.items():
        fig, axes = plt.subplots(1, 4, figsize=(15, 4), constrained_layout=True)
        axes[0].imshow(view_rgb(a)); axes[0].set_title('RGB display')
        for channel_index, ch_name in enumerate(('R', 'G', 'B')):
            axes[channel_index+1].imshow(view_channel(a, channel_index), cmap='gray', vmin=0, vmax=1)
            axes[channel_index+1].set_title(f'{ch_name}: within-image display')
        for ax in axes: ax.axis('off')
        fig.suptitle(f'{group}/{name} | channel identity NOT experimentally verified')
        safe_name = name.replace('.tif','').replace('+','_plus_').replace(' ','_')
        fig.savefig(out / f'channels_{group}_{safe_name}.png', dpi=140)
        plt.close(fig)

def main():
    parser = argparse.ArgumentParser(description='Build QC and RGB channel previews for three cell lines; no biological quantification.')
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    cache, rows = {}, []
    for label, group, names in CELL_FILES:
        for name in names:
            p = args.root / group / name
            record = {'cell_line_assumed':label, 'folder':group, 'file':name, 'exists':p.is_file(), 'status':'',
                      'height':'','width':'','channels':'','dtype':'','photometric':'','R_mean':'','G_mean':'','B_mean':'',
                      'interpretation':'RGB export; channel assignments unverified; not absolute HER2 mRNA expression'}
            if p.is_file():
                try:
                    a, photometric = read_tiff(p)
                    cache[(group,name)] = (a,photometric)
                    record.update({'status':'read_ok','height':a.shape[0],'width':a.shape[1], 'channels':a.shape[2],
                                   'dtype':str(a.dtype), 'photometric':photometric, 
                                   'R_mean':f'{np.mean(a[...,0]):.3f}', 'G_mean':f'{np.mean(a[...,1]):.3f}',
                                   'B_mean':f'{np.mean(a[...,2]):.3f}'})
                except Exception as e:
                    record['status'] = f'error: {e}'
            else:
                record['status']='missing'
            rows.append(record)
    with (args.out/'image_inventory.csv').open('w', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader(); writer.writerows(rows)
    plot_overview(cache, args.out)
    plot_channels(cache, args.out)
    print('Completed:',len(cache),'readable images of 9 expected')
    print('Saved:',args.out)
    print('Start by checking: three_cell_line_overview.png and channels_231_HER2-eu.png')
    print('WARNING: No mRNA spot quantification was performed; these are display RGB TIFF exports.')

if __name__ == '__main__':
    main()
