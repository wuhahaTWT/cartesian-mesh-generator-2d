#!/usr/bin/env python3
"""Physical free-slip reference, restart and invalid-input checks through the CLI."""
import argparse
from pathlib import Path
import sys
import tempfile
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools/verification'))
import verify_symmetry_flow as symmetry
p=argparse.ArgumentParser();p.add_argument('--cli',required=True);a=p.parse_args()
with tempfile.TemporaryDirectory(prefix='cm2d-symmetry-') as name:
    for scheme in ('upwind','limited-linear','face-limited-linear'):
        assert symmetry.study(Path(name)/scheme,a.cli,(4,),scheme)['valid']
