"""Run shearlib's built-in invariant checks (python lab/shearlib.py) under pytest.

One test for the whole selftest; the per-property tests are in test_units.py.
"""
import shearlib


def test_selftest():
    shearlib.selftest()
