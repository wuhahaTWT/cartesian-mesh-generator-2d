#!/usr/bin/env python3
"""Real synchronized passive-thermal cases, restart equivalence and native readback.

Each scalar step uses the newly converged flow flux. This is one-way constant-
property transport, not buoyancy or conjugate heat transfer. Serial 90 s/run cap.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import time
import verify_native_flow as native
import verify_scalar_transport as scalar
import verify_transient_flow as transient


def generate(root, mesh_cli, flow_cli, transport_cli):
    scalar.require(not root.exists(), 'output must be a fresh directory')
    root.mkdir(parents=True)
    summary=dict(valid=False,runs=[],spatial=[],temporal=[],circle=None,restart=None,carrierAudit=None,
                 scope='synchronized one-way passive thermal transport; fixed properties and thermal BC/source; no desktop or external CFD qualification',
                 executables={str(p):native.sha256_file(p) for p in (mesh_cli,flow_cli,transport_cli)})
    def execute(cmd,label):
        start=time.monotonic();record=native.run(cmd,root/'logs'/label,90.)
        record['wallSeconds']=time.monotonic()-start
        summary['runs'].append(record)
        scalar.require(record['returncode']==0 and not record['timedOut'],f'{label} failed; see recorded log')
    def vortex(mesh,out,dt,steps,diffusivity=.02,restart=None):
        cmd=[str(transport_cli),'--mesh',str(mesh),'--output',str(out),'--verification','thermal-vortex',
             '--dt',str(dt),'--steps',str(steps),'--flow-nu','.1','--diffusivity',str(diffusivity),
             '--flow-convection','limited-linear','--convection','limited-linear','--pressure-preconditioner','aggregation']
        if restart:cmd+=['--restart',str(restart)]
        execute(cmd,out.name)
        audit=scalar.verify(out)
        audit.update(prefix=str(out),dt=dt,time=json.loads(Path(str(out)+'.json').read_text())['time'])
        return audit
    try:
        meshes=[]
        for level in (4,5,6):
            prefix=root/'meshes'/f'l{level}'/'mesh';prefix.parent.mkdir(parents=True)
            request=native.Request(f'l{level}','manufactured',level,1/((1<<level)-2),.1,1.)
            execute(native.mesh_command(mesh_cli,prefix,request,root),f'mesh-l{level}')
            meshes.append(Path(str(prefix)+'.solver.cm2d'))
            summary['spatial'].append(vortex(meshes[-1],root/'spatial'/f'l{level}',.0025,20))
        for count in (5,10,20):
            summary['temporal'].append(vortex(meshes[1],root/'time'/f'n{count}',.1/count,count,.2))
        for key in ('spatial','temporal'):
            for coarse,fine in zip(summary[key],summary[key][1:]):
                ratio=coarse['dt']/fine['dt'] if key=='temporal' else math.sqrt(fine['cells']/coarse['cells'])
                fine['observedOrder']=math.log(coarse['l2Error']/fine['l2Error'])/math.log(ratio)
                scalar.require(fine['l2Error']<coarse['l2Error'],key+' error did not decrease')
        continuous=root/'restart'/'continuous';split=root/'restart'/'split';resumed=root/'restart'/'resumed'
        vortex(meshes[0],continuous,.01,4,.2)
        vortex(meshes[0],split,.01,2,.2)
        vortex(meshes[0],resumed,.01,2,.2,Path(str(split)+'.thermal.checkpoint'))
        matches={suffix:Path(str(continuous)+suffix).read_bytes()==Path(str(resumed)+suffix).read_bytes()
                 for suffix in ('.thermal.checkpoint','.carrier.checkpoint','.cells.csv','.faces.csv')}
        scalar.require(all(matches.values()),'restart differs from continuous calculation')
        summary['restart']=dict(byteIdentical=matches,continuousPrefix=str(continuous),resumedPrefix=str(resumed))
        flow_prefix=root/'flow-independent'/'vortex'
        execute([str(flow_cli),'--mesh',str(meshes[0]),'--output',str(flow_prefix),'--case','taylor-green',
                 '--time-step','.01','--steps','4','--nu','.1','--speed','1','--tolerance','1e-8',
                 '--convection','limited-linear','--pressure-preconditioner','aggregation'],'flow-independent')
        carrier_match=Path(str(flow_prefix)+'.checkpoint').read_bytes()==Path(str(continuous)+'.carrier.checkpoint').read_bytes()
        scalar.require(carrier_match,'thermal carrier differs from independently executed flow')
        summary['carrierAudit']=transient.verify(meshes[0],flow_prefix,root/'carrier-audit.json')
        scalar.require(summary['carrierAudit']['valid'],'independent flow equation audit failed')
        summary['carrierAudit']['bitwiseMatchesCoupledCarrier']=carrier_match
        prefix=root/'meshes'/'circle'/'mesh';prefix.parent.mkdir(parents=True)
        execute(native.mesh_command(mesh_cli,prefix,native.Request('circle','external',None,None,.1,1),root),'mesh-circle')
        mesh_path=Path(str(prefix)+'.solver.cm2d');mesh=native.read_cm2d(mesh_path)
        boundary=root/'circle-boundary.csv'
        with boundary.open('w') as f:
            writer=csv.writer(f,lineterminator='\n');writer.writerow(['face','type','value','inflowValue'])
            for e in mesh.edges:
                if e.neighbour<0:writer.writerow([e.id,'value' if e.patch==1 else 'flux',1 if e.patch==1 else 0,'' if e.patch==1 else 0])
        out=root/'circle'/'result'
        execute([str(transport_cli),'--mesh',str(mesh_path),'--output',str(out),'--evolve-flow','external','--boundary',str(boundary),
                 '--flow-nu','.1','--flow-speed','1','--flow-convection','limited-linear','--pressure-preconditioner','aggregation',
                 '--diffusivity','.1','--initial','0','--dt','.05','--steps','4'],'heated-circle')
        summary['circle']=scalar.verify(out);summary['circle']['prefix']=str(out)
        summary['valid']=True
    except (ValueError,OSError,KeyError) as exc:
        summary['error']=str(exc)
    (root/'audit.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    return summary

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate',type=Path,required=True)
    parser.add_argument('--mesh-cli',type=Path,default=Path('build/cartmesh2d_cli'))
    parser.add_argument('--flow-cli',type=Path,default=Path('build/cartmesh2d_flow_cli'))
    parser.add_argument('--transport-cli',type=Path,default=Path('build/cartmesh2d_transport_cli'))
    args=parser.parse_args()
    report=generate(args.generate.resolve(),args.mesh_cli.resolve(),args.flow_cli.resolve(),args.transport_cli.resolve())
    print(json.dumps(dict(valid=report['valid'],error=report.get('error'),runs=len(report['runs']))))
    raise SystemExit(0 if report['valid'] else 1)
