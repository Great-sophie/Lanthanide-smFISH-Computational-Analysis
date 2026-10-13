#!/usr/bin/env python3
"""Audit TrackMate spot exports per SKBR3 cell and optionally compare matched ImageJ counts.

Each gradient-number CSV is one cell, per experimenter confirmation.
All spot counts remain algorithm-specific candidate detections.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

GRADIENTS = ('25%', '50%', '75%', '100%')
STEM = re.compile(r'^(25%|50%|75%|100%)-(\d+)$')
QUALITY_THRESHOLDS = (0.5, 1.0, 2.0)


def read_trackmate(path: Path) -> dict:
    with path.open(encoding='utf-8-sig', errors='replace', newline='') as f:
        rows = list(csv.reader(f))
    if len(rows) < 4:
        raise ValueError('fewer than four TrackMate header/metadata rows')
    head = [str(s).strip().upper() for s in rows[0]]
    expected = ('ID', 'QUALITY', 'POSITION_X', 'POSITION_Y')
    if not set(expected).issubset(head):
        raise ValueError('missing TrackMate ID/QUALITY/X/Y columns')
    ix = {field: head.index(field) for field in expected}
    ids, qs = [], []
    footer_rows = malformed = 0
    unit_row = rows[3]
    units = {k: unit_row[ix[k]].strip() if len(unit_row) > ix[k] else '' for k in ('POSITION_X', 'POSITION_Y')}
    for row in rows[4:]:
        get = lambda field: row[ix[field]].strip() if len(row) > ix[field] else ''
        sid, qraw, xraw, yraw = [get(k) for k in expected]
        # TrackMate sometimes places summary statistics in columns other than
        # ID/QUALITY/X/Y at the end: they are not malformed Spot detections.
        if not any((sid, qraw, xraw, yraw)):
            footer_rows += 1
            continue
        try:
            q, x, y = (float(qraw), float(xraw), float(yraw))
            if not sid or not all(math.isfinite(v) for v in (q, x, y)):
                raise ValueError('incomplete/nonfinite spot')
        except (ValueError, TypeError):
            malformed += 1
            continue
        ids.append(sid)
        qs.append(q)
    if len(ids) != len(set(ids)):
        raise ValueError('duplicate Spot IDs within CSV')
    return dict(spot_count=len(ids), qualities=np.array(qs, dtype=float),
                footer_rows=footer_rows, malformed_rows=malformed,
                coord_units=f"X={units['POSITION_X']}; Y={units['POSITION_Y']}")


def audit(root: Path) -> list[dict]:
    records = []
    for group in GRADIENTS:
        folder = root / group
        if not folder.is_dir():
            raise FileNotFoundError(f'Missing gradient directory: {folder}')
        for path in sorted(folder.glob('*.csv')):
            match = STEM.fullmatch(path.stem)
            fname_group = match.group(1) if match else ''
            rec = dict(gradient=group, filename_gradient=fname_group, cell_id=path.stem,
                       csv_path=str(path.relative_to(root)), status='excluded', reason='',
                       spot_count='', median_quality='', footer_rows='', malformed_rows='',
                       coord_units='', has_matching_tiff=path.with_suffix('.tif').is_file(),
                       sha256=hashlib.sha256(path.read_bytes()).hexdigest())
            for q in QUALITY_THRESHOLDS:
                k=str(q).replace('.','p')
                rec[f'count_quality_ge_{k}']=''
                rec[f'fraction_quality_ge_{k}']=''
            if not match:
                rec['reason'] = 'unrecognized filename; manual identity review'
            elif fname_group != group:
                rec['reason'] = 'filename gradient conflicts with directory; not automatically reassigned'
            else:
                try:
                    info = read_trackmate(path)
                    for key in ('spot_count', 'median_quality', 'footer_rows', 'malformed_rows', 'coord_units'):
                        if key == 'median_quality':
                            rec[key] = float(np.median(info['qualities'])) if info['spot_count'] else ''
                        else:
                            rec[key] = info[key]
                    if info['malformed_rows']:
                        rec['reason'] = f"{info['malformed_rows']} malformed spot rows; manual review required"
                    elif not info['spot_count']:
                        rec['reason'] = 'no valid spots; zero vs export issue unverified'
                    else:
                        rec['status'] = 'included'
                        if not rec['has_matching_tiff']:
                            rec['reason'] = 'no same-stem TIFF; count still valid as confirmed single-cell CSV'
                    for q in QUALITY_THRESHOLDS:
                        k=str(q).replace('.','p')
                        n=int(np.count_nonzero(info['qualities'] >= q))
                        rec[f'count_quality_ge_{k}']=n
                        rec[f'fraction_quality_ge_{k}']=n/info['spot_count'] if info['spot_count'] else ''
                except Exception as exc:
                    rec['reason'] = f'parse error: {exc}'
            records.append(rec)
    seen = {}
    for r in records:
        if r['status'] != 'included':
            continue
        if r['sha256'] in seen:
            r['status'] = 'excluded'
            r['reason'] = f'exact duplicate bytes of {seen[r["sha256"]]}; investigate'
        else:
            seen[r['sha256']] = r['csv_path']
    return records


def summaries(records):
    summaries=[]
    for group in GRADIENTS:
        included=[r for r in records if r['gradient']==group and r['status']=='included']
        values=np.array([r['spot_count'] for r in included], dtype=float)
        out=dict(gradient=group, included_cells=len(included),
                 excluded_files=sum(r['gradient']==group and r['status']!='included' for r in records),
                 mean_spots=float(np.mean(values)) if len(values) else '',
                 median_spots=float(np.median(values)) if len(values) else '',
                 sd_spots=float(np.std(values,ddof=1)) if len(values)>1 else '',
                 q1_spots=float(np.percentile(values,25)) if len(values) else '',
                 q3_spots=float(np.percentile(values,75)) if len(values) else '')
        for q in QUALITY_THRESHOLDS:
            k=str(q).replace('.','p')
            v=[r[f'fraction_quality_ge_{k}'] for r in included]
            out[f'median_fraction_quality_ge_{k}']=float(np.median(v)) if v else ''
        summaries.append(out)
    return summaries


def write_csv(path:Path, records:list[dict]) -> None:
    if not records:
        return
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        wr=csv.DictWriter(f, fieldnames=list(records[0]));wr.writeheader();wr.writerows(records)


def plot_trackmate(records, outfile):
    fig,ax=plt.subplots(figsize=(9.5,5.7))
    for i,g in enumerate(GRADIENTS):
        r=sorted((r for r in records if r['gradient']==g and r['status']=='included'),key=lambda r:r['cell_id'])
        y=np.array([r['spot_count'] for r in r],float)
        if len(y):
            offsets=np.linspace(-.16,.16,len(y)) if len(y)>1 else [0]
            ax.scatter(i+np.array(offsets),y,s=45,alpha=.85)
            med=float(np.median(y))
            ax.hlines(med,i-.23,i+.23,color='black',lw=2.6)
            ax.annotate(f'n={len(y)}; median={med:g}',(i,max(y)),xytext=(0,12),textcoords='offset points',ha='center',fontsize=9)
    ax.set_xticks(range(4),GRADIENTS)
    ax.set_xlim(-.55,3.55)
    ax.set_xlabel('dUTP-N3 labeling fraction')
    ax.set_ylabel('TrackMate candidate spots per cell')
    ax.set_title('SKBR3 HER2 smFISH — historical TrackMate counts')
    ax.grid(axis='y',alpha=.25)
    fig.text(.5,.018,'Single-cell measurement ≠ independent biological replicate. TrackMate counts are not confirmed individual HER2 transcripts.',ha='center',fontsize=8.5)
    fig.tight_layout(rect=(0,.055,1,1));fig.savefig(outfile,dpi=190);plt.close(fig)


def plot_quality(records,outfile):
    fig,axes=plt.subplots(1,2,figsize=(11,4.5))
    for q in QUALITY_THRESHOLDS:
        k=str(q).replace('.','p')
        med_counts=[];med_fracs=[]
        for g in GRADIENTS:
            rs=[r for r in records if r['status']=='included' and r['gradient']==g]
            med_counts.append(float(np.median([r[f'count_quality_ge_{k}'] for r in rs])) if rs else np.nan)
            med_fracs.append(100*float(np.median([r[f'fraction_quality_ge_{k}'] for r in rs])) if rs else np.nan)
        axes[0].plot(range(4),med_counts,marker='o',label=f'Q ≥ {q:g}')
        axes[1].plot(range(4),med_fracs,marker='o',label=f'Q ≥ {q:g}')
    for ax in axes:
        ax.set_xticks(range(4),GRADIENTS);ax.set_xlabel('dUTP-N3 labeling fraction')
        ax.grid(axis='y',alpha=.25);ax.legend()
    axes[0].set_ylabel('Median retained candidates / cell')
    axes[1].set_ylabel('Median retained fraction (%)')
    axes[1].set_ylim(0,105)
    fig.suptitle('TrackMate detector QUALITY threshold sensitivity')
    fig.text(.5,.015,'QUALITY is an uncalibrated detector score; compare across runs only if settings are consistent.',ha='center',fontsize=8.5)
    fig.tight_layout(rect=(0,.06,1,.95));fig.savefig(outfile,dpi=180);plt.close(fig)


def read_imagej(path:Path)->list[dict]:
    with path.open(encoding='utf-8-sig',newline='') as f:
        rows=list(csv.DictReader(f))
    for expected in ('gradient','cell_id','imagej_spots'):
        if not rows or expected not in rows[0]:
            raise ValueError('ImageJ file needs columns: gradient,cell_id,imagej_spots')
    out=[];seen=set()
    for i,r in enumerate(rows,2):
        g=r['gradient'].strip(); c=r['cell_id'].strip()
        if g not in GRADIENTS or not STEM.fullmatch(c) or not c.startswith(g+'-'):
            raise ValueError(f'ImageJ row {i}: invalid gradient/cell_id')
        try:
            n=float(r['imagej_spots'])
        except (ValueError,TypeError):
            raise ValueError(f'ImageJ row {i}: nonnumeric spot count')
        if not math.isfinite(n) or n<0 or not n.is_integer():
            raise ValueError(f'ImageJ row {i}: spot count must be a nonnegative integer')
        if (g,c) in seen:
            raise ValueError(f'Duplicate ImageJ cell_id: {c}')
        seen.add((g,c))
        out.append(dict(gradient=g,cell_id=c,imagej_spots=int(n)))
    return out


def matched_methods(trackmate, imagej):
    by_cell={(r['gradient'],r['cell_id']):r for r in trackmate if r['status']=='included'}
    records=[]
    for r in imagej:
        key=(r['gradient'],r['cell_id'])
        if key not in by_cell:
            continue
        tr=int(by_cell[key]['spot_count']);ij=r['imagej_spots']
        records.append(dict(gradient=r['gradient'],cell_id=r['cell_id'],imagej_spots=ij,
                            trackmate_spots=tr,trackmate_minus_imagej=tr-ij,
                            mean_of_methods=(tr+ij)/2))
    return records


def plot_matched(pairs,outfile):
    fig,axs=plt.subplots(1,2,figsize=(10,4.7))
    xs=np.array([r['imagej_spots'] for r in pairs],dtype=float)
    ys=np.array([r['trackmate_spots'] for r in pairs],dtype=float)
    if len(xs):
        for g in GRADIENTS:
            group=[r for r in pairs if r['gradient']==g]
            if group:
                axs[0].scatter([r['imagej_spots'] for r in group], [r['trackmate_spots'] for r in group],label=g,s=45)
                axs[1].scatter([r['mean_of_methods'] for r in group], [r['trackmate_minus_imagej'] for r in group],s=45,label=g)
        upper=max(xs.max(),ys.max())*1.08+1
        axs[0].plot([0,upper],[0,upper],'--',color='gray',lw=1)
        axs[0].set_xlim(0,upper);axs[0].set_ylim(0,upper)
        axs[1].axhline(0,linestyle='--',color='gray',lw=1)
    axs[0].set(xlabel='ImageJ candidates per same cell',ylabel='TrackMate candidates per same cell',title='Matched cell counts')
    axs[1].set(xlabel='Mean of methods',ylabel='TrackMate − ImageJ',title='Difference vs mean (descriptive)')
    for a in axs:a.grid(alpha=.2);a.legend(loc='best')
    fig.suptitle(f'Algorithm comparison: n={len(pairs)} confirmed matched cells')
    fig.tight_layout();fig.savefig(outfile,dpi=180);plt.close(fig)


def main(argv=None):
    p=argparse.ArgumentParser(description='SKBR3 single-cell TrackMate QC and optional paired ImageJ comparison')
    p.add_argument('--root',required=True,type=Path,help='fig3 0710 containing 25%%, 50%%, 75%%, 100%% subfolders')
    p.add_argument('--out',required=True,type=Path)
    p.add_argument('--imagej',type=Path,help='OPTIONAL CSV: gradient,cell_id,imagej_spots; only if SAME cells confirmed')
    args=p.parse_args(argv)
    args.out.mkdir(parents=True,exist_ok=True)
    records=audit(args.root)
    write_csv(args.out/'trackmate_cell_audit.csv',records)
    stats=summaries(records)
    write_csv(args.out/'trackmate_gradient_summary.csv',stats)
    write_csv(args.out/'excluded_and_review.csv',[r for r in records if r['status']!='included' or r['malformed_rows'] not in ('',0) or 'no same-stem' in r['reason']])
    plot_trackmate(records,args.out/'trackmate_spots_per_cell.png')
    plot_quality(records,args.out/'trackmate_quality_sensitivity.png')
    included=sum(r['status']=='included' for r in records)
    notes=[f'TrackMate CSV files audited: {len(records)}; included {included}; excluded {len(records)-included}.',
           'Blank/metadata footer rows excluded from candidate-spot counts (not considered malformed).',
           'Folder/filename gradient conflicts remain excluded; no automatic reassignment.',
           'TrackMate counts are exploratory; DO NOT claim ImageJ-vs-TrackMate algorithm comparison without matched ImageJ counts.',
           'TrackMate QUALITY is not calibrated SNR; detector settings and independent biological replicates must be documented.',
           'Groups may contain unequal numbers of cells; no significance testing is performed.',
           '', 'GRADIENT SUMMARY:']
    for s in stats:
        notes.append(f"{s['gradient']}: n={s['included_cells']}, mean={s['mean_spots']}, median={s['median_spots']}")
    if args.imagej:
        imagej=read_imagej(args.imagej)
        pairs=matched_methods(records,imagej)
        write_csv(args.out/'matched_imagej_trackmate.csv',pairs)
        notes.extend(['',f'ImageJ rows supplied: {len(imagej)}; confirmed same-cell IDs matched to TrackMate: {len(pairs)}.',
                      'Matched results are descriptive; verify same raw images, ROIs, preprocessing and reference conditions.'])
        if pairs:
            plot_matched(pairs,args.out/'matched_method_comparison.png')
        else:
            notes.append('No matches; method comparison not generated.')
    (args.out/'interpretation.txt').write_text('\n'.join(notes)+'\n',encoding='utf-8')
    print('\n'.join(notes))
    print('Results:',args.out)
    return records,stats

if __name__=='__main__':
    main()
