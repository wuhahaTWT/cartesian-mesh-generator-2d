#!/usr/bin/env python3
"""Rotating smooth-wall accuracy, template rejection, covariance and restart."""
import argparse
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools/verification'))
import verify_rotating_annulus as study
import verify_native_flow as native
import verify_transient_flow as transient
from verify_named_channel import transform_mesh

def run(args,success=True):
    r=subprocess.run(list(map(str,args)),text=True,capture_output=True,timeout=30)
    assert (r.returncode==0)==success,r.stdout+r.stderr
    return r

def main(cli):
    with tempfile.TemporaryDirectory(prefix='cm2d-rotating-') as tmp:
        root=Path(tmp);mesh=root/'ring.solver.cm2d';study.write_polar_mesh(mesh,8,64)
        bc=root/'input.boundaries'
        def template(m,path,success=True):return run([cli,'--mesh',m,'--case','annulus','--speed',.5,'--export-boundaries',path],success)
        template(mesh,bc)
        assert bc.read_text().count(' smooth-moving-wall ')==64
        assert bc.read_text().count(' wall ')==64
        def solve(prefix,m=mesh,boundary=bc,extra=(),success=True):
            return run([cli,'--mesh',m,'--case','custom','--boundary',boundary,'--output',prefix,
                '--nu',.1,'--speed',.5,'--convection','face-limited-linear','--pressure-preconditioner','aggregation',
                '--max-iterations',2000,'--tolerance',1e-9,*extra],success)
        def audit(prefix,m=mesh):
            a=native.verify_case(m,prefix,'custom',.1,.5,native.argument_parser().parse_args(['--max-iterations','2000']))
            assert a['valid'],a['issues']
        base=root/'base';solve(base);audit(base);metrics=study.metrics(mesh,base)
        assert metrics['velocityL2']<.01 and metrics['pressureL2']<.005,metrics
        assert max(metrics['rotorTorqueError'],metrics['housingTorqueError'])<.02,metrics
        # The legacy explicit constant trace remains distinct and restart-incompatible.
        old=root/'constant.boundaries';old.write_text(bc.read_text().replace('smooth-moving-wall','moving-wall'))
        constant=root/'constant';solve(constant,boundary=old);audit(constant)
        assert metrics['velocityL2']<.4*study.metrics(mesh,constant)['velocityL2']
        turned=root/'rotated.solver.cm2d';transform_mesh(mesh,turned,.63)
        turned_bc=root/'turned.boundaries';template(turned,turned_bc)
        rotated=root/'turned';solve(rotated,m=turned,boundary=turned_bc);audit(rotated,turned)
        def cells(prefix):return native.load_csv(Path(str(prefix)+'.cells.csv'),('u','v','p'))
        c,s=math.cos(.63),math.sin(.63)
        for a,b in zip(cells(base),cells(rotated)):
            u,v=float(a['u']),float(a['v'])
            assert math.hypot(c*u-s*v-float(b['u']),s*u+c*v-float(b['v']))<1e-7
            assert abs(float(a['p'])-float(b['p']))<1e-7
        whole,first,resumed=(root/k for k in ('whole','first','resumed'))
        solve(whole,extra=['--time-step',.02,'--steps',2]);solve(first,extra=['--time-step',.02,'--steps',1])
        solve(resumed,extra=['--time-step',.02,'--steps',1,'--restart',str(first)+'.checkpoint'])
        assert Path(str(whole)+'.checkpoint').read_bytes()==Path(str(resumed)+'.checkpoint').read_bytes()
        assert ' smooth-moving-wall ' in Path(str(whole)+'.checkpoint').read_text()
        for prefix in (whole,resumed):assert transient.verify(mesh,prefix,root/(prefix.name+'-audit.json'))['valid']
        rejection=solve(root/'wrong-kind',boundary=old,extra=['--time-step',.02,'--steps',1,'--restart',str(first)+'.checkpoint'],success=False)
        assert 'boundary' in rejection.stderr.lower(),rejection.stderr
        # Squares and odd polygons are explicitly outside this template's scope.
        for sides in (4,31):
            rejected_mesh=root/f'rejected-{sides}.solver.cm2d';study.write_polar_mesh(rejected_mesh,4,sides)
            rejected=template(rejected_mesh,root/'rejected.boundaries',False)
            assert 'annulus' in rejected.stderr.lower(),rejected.stderr
        # A valid deformed ring with eccentric walls must not be misclassified.
        eccentric=root/'eccentric.solver.cm2d';lines=mesh.read_text().splitlines();nv=int(lines[1].split()[1]);vertices=[]
        for i in range(2,nv+2):
            ident,x,y=lines[i].split();x,y=float(x),float(y);r=math.hypot(x,y);x+=.08*(1-r)
            vertices.append((x,y));lines[i]=f'{ident} {x:.17g} {y:.17g}'
        cell_start=next(i for i,l in enumerate(lines) if l.startswith('CELLS '))
        for i in range(cell_start+1,cell_start+1+int(lines[cell_start].split()[1])):
            t=lines[i].split();ids=list(map(int,t[5:5+int(t[4])]))
            t[3]=repr(.5*math.fsum(vertices[a][0]*vertices[b][1]-vertices[b][0]*vertices[a][1] for a,b in zip(ids,ids[1:]+ids[:1])))
            lines[i]=' '.join(t)
        eccentric.write_text('\n'.join(lines)+'\n')
        rejected=template(eccentric,root/'eccentric.boundaries',False)
        assert 'concentric' in rejected.stderr.lower(),rejected.stderr
        print(json.dumps(dict(valid=True,metrics=metrics,restartIdentical=True,rotation=.63)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cli',type=Path,required=True);a=p.parse_args();main(a.cli.resolve())
