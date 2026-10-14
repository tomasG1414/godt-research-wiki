"""Run shearlib's invariant checks under pytest (so CI reports them by name)."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "lab"))

import matplotlib

matplotlib.use("Agg")

import shearlib  # noqa: E402


def test_selftest():
    shearlib.selftest()
