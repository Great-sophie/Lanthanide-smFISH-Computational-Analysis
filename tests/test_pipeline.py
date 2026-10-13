import csv
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "imaging_maps"))
sys.path.insert(0, str(ROOT / "spot_counts"))

from chapter3_mrna_heatmap_global_only import local_signal, tile_mean, robust_high
from chapter3_spot_scatter_stats import read_spots
from method_robustness import read_trackmate

HEADER = "LABEL,ID,TRACK_ID,QUALITY,POSITION_X,POSITION_Y,POSITION_Z,POSITION_T,FRAME,RADIUS,VISIBILITY,MANUAL_SPOT_COLOR\n"
DESC1 = "Label,Spot ID,Track ID,Quality,X,Y,Z,T,Frame,Radius,Visibility,Manual spot color\n"
DESC2 = "Label,Spot ID,Track ID,Quality,X,Y,Z,T,Frame,R,Visibility,Spot color\n"
UNITS = ",,,(quality),(pixel),(pixel),(pixel),(sec),,(pixel),,\n"


def test_local_contrast_nonnegative():
    x = np.ones((80, 80), dtype=np.uint8) * 5
    x[35:42, 35:42] = 80
    y = local_signal(x, 2, 20)
    assert y.shape == x.shape
    assert y.min() >= 0
    assert y[38, 38] > y[0, 0]
    assert robust_high(y, 99.5) > 0


def test_tile_mean_constant():
    x = np.ones((11, 10), dtype=np.float32) * .25
    out = tile_mean(x, 4)
    assert out.shape == (3, 3)
    assert np.allclose(out, .25)


def test_trackmate_footer_ignored(tmp_path):
    p = tmp_path / '25%-1.csv'
    p.write_text(HEADER + DESC1 + DESC2 + UNITS
                 + 'ID1,1,,0.9,10,10,0,0,0,5,1,\n'
                 + 'ID2,2,,1.4,15,14,0,0,0,5,1,\n'
                 + ',,,,,,,,,,,,,summary,22\n', encoding='utf-8')
    count, blanks, invalid = read_spots(p)
    assert count == 2 and blanks == 1 and invalid == 0
    info = read_trackmate(p)
    assert info['spot_count'] == 2 and info['footer_rows'] == 1
