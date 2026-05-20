#!/usr/bin/env python3
"""Compatibility wrapper for running EasyScan_HEP from the source tree."""

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
for path in (ROOT, ROOT / "src"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from easyscan_hep.cli import main


if __name__ == "__main__":
    main()
