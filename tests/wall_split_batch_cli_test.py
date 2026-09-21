#!/usr/bin/env python3
"""Actual nozzle wall defects: full gates and bounded global rebuild work."""
import argparse,json,subprocess,sys,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools/verification'))
import verify_native_flow as reader
p=argparse.ArgumentParser();p.add_argument('--cli',type=Path,required=True);args=p.parse_args()
with tempfile.TemporaryDirectory(prefix='cm2d-wall-batch-') as directory:
    root=Path(directory);prefix=root/'nozzle';source=ROOT/'examples/complex/nozzle_profile.xy'
    result=subprocess.run([str(args.cli.resolve()),str(source),str(prefix),'7','0.03333333333333333','.1','interior',str(root/'foam'),'7','0'],capture_output=True,text=True,timeout=60)
    assert result.returncode==0,result.stdout+result.stderr
    profile=dict(line.split('=',1) for line in result.stdout.splitlines() if '=' in line)
    assert profile['solver_quality']=='PASS'
    # Counts bound algorithmic work, not machine-dependent wall time. The
    # missing single-cell proposals needed 1760 full rebuilds on this fixture.
    assert int(profile['solver_profile_global_topology_rebuild_calls'])<500
    assert int(profile['solver_profile_repartition_iterations'])<4
    mesh=reader.read_cm2d(prefix.with_suffix('.solver.cm2d'));measured=reader.measure(mesh,1e-11,1e-9)
    assert not measured.issues,measured.issues
    points=[tuple(map(float,line.split())) for line in source.read_text().splitlines() if line.strip()]
    expected=abs(sum(a[0]*b[1]-a[1]*b[0] for a,b in zip(points,points[1:]+points[:1])))/2
    assert abs(measured.total_area-expected)<1e-10
    # Every exported physical face remains on an original input segment.
    def on_segment(p,a,b):
        dx,dy=b[0]-a[0],b[1]-a[1];cross=(p[0]-a[0])*dy-(p[1]-a[1])*dx
        return abs(cross)<1e-10 and min(a[0],b[0])-1e-10<=p[0]<=max(a[0],b[0])+1e-10 and min(a[1],b[1])-1e-10<=p[1]<=max(a[1],b[1])+1e-10
    boundary_count=0
    for edge in mesh.edges:
        if edge.neighbour>=0:continue
        boundary_count+=1
        a,b=mesh.vertices[edge.v0],mesh.vertices[edge.v1]
        assert any(on_segment(a,x,y) and on_segment(b,x,y) for x,y in zip(points,points[1:]+points[:1]))
    assert boundary_count>0
print('Actual nozzle wall batch: independent topology/area/original boundary and full solver gates pass')
