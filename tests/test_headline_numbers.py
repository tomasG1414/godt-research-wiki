"""The numbers README.md quotes must be the numbers the scripts print.

Runs the fast experiments exactly as a reader would (python expNN.py in lab/)
and checks the printed values behind the README's headline claims. ~50 s.
If a change moves one of these numbers, update the wiki page and this test
together. Skip locally with: pytest -m "not slow"
"""
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

LAB = Path(__file__).resolve().parent.parent / "lab"
pytestmark = pytest.mark.slow


def run(script):
    env = dict(os.environ, MPLBACKEND="Agg")
    res = subprocess.run([sys.executable, script], cwd=LAB, env=env,
                         capture_output=True, text=True, timeout=900)
    assert res.returncode == 0, res.stdout[-2000:] + res.stderr[-2000:]
    return res.stdout


def row(text, label):
    """The rest of the first output line containing `label`."""
    for line in text.splitlines():
        if label in line:
            return line.split(label, 1)[1]
    raise AssertionError(f"{label!r} not found in the script output")


def floats(s):
    return [float(x) for x in re.findall(r"-?\d+\.\d+", s)]


def after(s, key):
    m = re.search(re.escape(key) + r"\s*(-?\d+\.\d+)", s)
    assert m, f"{key!r} not found in {s!r}"
    return float(m.group(1))


@pytest.fixture(scope="module")
def exp01():
    return run("exp01_shear_nulls.py")


@pytest.fixture(scope="module")
def exp02():
    return run("exp02_shear_recovery.py")


@pytest.fixture(scope="module")
def exp03():
    return run("exp03_ewald_fisher.py")


@pytest.fixture(scope="module")
def exp05b():
    return run("exp05b_weighted_average.py")


def test_exp01_published_shear_has_three_in_band_nulls(exp01):
    assert "band |K_x| <= 31.19 rad/um" in exp01
    assert "7.85, 15.71, 23.56" in exp01


def test_exp02_band_lost_along_the_shear_axis_and_image_gain(exp02):
    # finite-shear-null-experiment.md sections 3 and 4; lab/README.md "3.56 %", "8.5x"
    assert floats(row(exp02, "d=0.8 (published exp)"))[2] == 3.56
    assert floats(row(exp02, "d=0.8+0.55 "))[2] == 0.00
    assert floats(row(exp02, "x+y at 0.8 and 0.55"))[:2] == [0.0180, 8.5]


def test_exp03_lost_volume_in_the_readme_figure(exp03):
    # README caption: 0.8 um loses 3 %, 0.1 um 2.5 %, {0.8, 0.55} x + y zero
    assert floats(row(exp03, "0.8 x  (published exp)"))[1:3] == [0.0342, 3.02]
    assert floats(row(exp03, "0.1 x  (GODT sim shear)"))[2] == 2.54
    assert floats(row(exp03, "{0.8, 0.55} x + {0.8, 0.55} y"))[2] == 0.00
    assert floats(row(exp03, "an m != 0 null plane :"))[0] == 83.6
    assert floats(row(exp03, "the blind plane K_x=0:"))[0] == 16.4
    assert after(row(exp03, "field  25.6 um"), "proj:") == 76.6


def test_exp05b_fusion_against_equal_weight_averaging(exp05b):
    # README "10-21x" (upper end); multishear-fusion-design-solution.md section 5
    r = row(exp05b, "{0.8,0.55}x+y")
    wiener, naive = after(r, "wiener"), after(r, "naive")
    assert (wiener, naive, after(r, "prec-weighted avg")) == (0.0078, 0.1607, 0.0079)
    assert round(naive / wiener) == 21
