#!/usr/bin/env python3
"""Reproduce the PD calibration tables from the recorded partitions.

Run: python calibrate.py
Requires the parent calibration directory; uses only the Python standard library.
"""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from refresh_saved_calibration import run


if __name__ == '__main__':
    run(('pd',))
