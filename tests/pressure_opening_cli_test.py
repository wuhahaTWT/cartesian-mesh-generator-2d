#!/usr/bin/env python3
"""Independent pressure-only startup, stationary equilibrium and restart audit."""
import argparse
from pathlib import Path
import sys
import tempfile
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools/verification'))
import verify_pressure_openings as opening
parser=argparse.ArgumentParser();parser.add_argument('--cli',required=True);args=parser.parse_args()
with tempfile.TemporaryDirectory(prefix='cm2d-pressure-opening-') as name:
    for scheme in ('upwind','limited-linear','face-limited-linear'):
        assert opening.study(Path(name)/scheme,args.cli,(8,),scheme)['valid']
