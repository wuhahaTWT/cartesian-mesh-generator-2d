#!/usr/bin/env python3
"""Independent field, conservation, analytic profile and unsupported-mode checks."""
import argparse,csv,json,subprocess,sys,tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools/verification'))
import verify_native_flow as native
from verify_pressure_openings import rectangle,pressure_boundaries
from verify_symmetry_flow import symmetry_file
parser=argparse.ArgumentParser();parser.add_argument('--cli',required=True);args=parser.parse_args()
with tempfile.TemporaryDirectory(prefix='cm2d-anderson-') as directory:
    root=Path(directory);mesh=root/'mesh.solver.cm2d';rectangle(mesh,4)
    template=root/'template';pressure=root/'pressure';bc=root/'boundaries'
    subprocess.run([args.cli,'--mesh',str(mesh),'--case','duct','--export-boundaries',str(template)],check=True,capture_output=True)
    pressure_boundaries(template,pressure,1.2,0);symmetry_file(pressure,bc)
    for scheme in ['upwind','limited-linear','face-limited-linear']:
        fields={};histories={}
        for mode in ['none','anderson']:
            prefix=root/(scheme+'-'+mode)
            command=[args.cli,'--mesh',str(mesh),'--boundary',str(bc),'--case','custom','--output',str(prefix),
                '--nu','.1','--speed','1','--tolerance','1e-11','--max-iterations','5000','--convection',scheme,
                '--pressure-preconditioner','aggregation','--steady-acceleration',mode]
            subprocess.run(command,check=True,capture_output=True,timeout=30)
            report=native.verify_case(mesh,prefix,'custom',.1,1,native.argument_parser().parse_args(['--max-iterations','5000']))
            assert report['valid'],report['issues']
            summary=json.loads(Path(str(prefix)+'.json').read_text());assert summary['converged'] and summary['steadyAcceleration']==mode
            assert summary['accelerationAccepted']+summary['accelerationRejected']==summary['accelerationCandidates']
            if mode=='anderson':assert summary['accelerationAccepted']>0
            else:assert summary['accelerationCandidates']==0
            fields[mode]=list(csv.DictReader(Path(str(prefix)+'.cells.csv').open()));histories[mode]=summary['iterations']
        for a,b in zip(fields['none'],fields['anderson']):
            assert max(abs(float(a[k])-float(b[k])) for k in ['u','v','p'])<2e-8
        assert histories['anderson']<histories['none'],histories
    for extra in [['--time-step','.02','--steps','1'],['--steady-acceleration','unknown']]:
        bad=subprocess.run(command+extra,capture_output=True,text=True,timeout=30)
        assert bad.returncode!=0 and 'steady-acceleration' in bad.stderr,bad
