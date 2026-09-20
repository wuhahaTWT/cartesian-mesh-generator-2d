#!/usr/bin/env python3
"""Plot actual external meshes/fields and integrated-traction sensitivity."""
import argparse
import json
import math
import sys
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'tools/verification'))
import verify_native_flow as native


def main(args):
    source=args.study/'study.json';report=json.loads(source.read_text())
    if report.get('valid') is not True or report.get('accuracyQualification')!='not-qualified' or len(report['cases'])!=3:
        raise ValueError('Expected three fully audited cases without an accuracy qualification claim')
    for name,digest in report['artifactSha256'].items():
        path=Path(name)
        if path.is_absolute() or '..' in path.parts or native.sha256_file(args.study/path)!=digest:
            raise ValueError('Source artifact digest differs: '+name)
    fields=[];meshes=[]
    for case in report['cases']:
        folder=args.study/case['label'];mesh_path=folder/'mesh.solver.cm2d'
        audit=native.verify_case(mesh_path,folder/'flow','external',.1,1,native.argument_parser().parse_args(['--max-iterations','4000']))
        if not audit['valid']:raise ValueError('Independent final-field audit failed')
        mesh=native.read_cm2d(mesh_path);measured=native.measure(mesh,1e-11,1e-9)
        fields.append(native.read_cells(folder/'flow.cells.csv',mesh,measured,'external'));meshes.append(mesh)
    fig,axes=plt.subplots(2,2,figsize=(14,9),layout='constrained');limit=max(math.hypot(c['u'],c['v']) for field in fields for c in field)
    for ax,mesh,field,case in zip(axes.flat,meshes,fields,report['cases']):
        polygons=[[mesh.vertices[v] for v in c.vertices] for c in mesh.cells]
        colour=PolyCollection(polygons,array=[math.hypot(c['u'],c['v']) for c in field],cmap='viridis',edgecolors='#5c666b',linewidths=.2)
        colour.set_clim(0,limit);ax.add_collection(colour)
        ax.set(xlim=(-2,7),ylim=(-3,3),aspect='equal',xlabel='x [m]',ylabel='y [m]',title=f"{case['cells']:,} actual Cut-cell cells; wall h/Lref={case['wallRelativeSize']:g}")
        fig.colorbar(colour,ax=ax,label='Speed [m/s]',shrink=.85)
    ax=axes[1,1];h=[c['wallRelativeSize'] for c in report['cases']];drag=[c['summary']['forceX'] for c in report['cases']]
    # Uref=1 and D=2 make 2*Fx/(Uref^2*D) numerically equal to Fx.
    ax.semilogx(h,drag,'o-',label='Stopping tolerance 1e-8')
    ax.scatter([h[-1]],[report['iterationSensitivity']['tightSummary']['forceX']],marker='x',s=60,label='Finest grid, tolerance 1e-10')
    for x,y,c in zip(h,drag,report['cases']):ax.annotate(f"{c['cells']:,} cells\nCd={y:.6f}",(x,y),xytext=(4,10),textcoords='offset points',fontsize=9)
    ax.margins(.3);ax.invert_xaxis();ax.grid(alpha=.25)
    ax.set(xlabel='Requested wall h/Lref (background refined together)',ylabel='Cd = 2 Fx / (Uref² D)',title='Integrated traction changes with resolution')
    ax.legend(fontsize=9)
    delta=report['iterationSensitivity'];ax.text(.03,.03,f"Finest tolerance change:\nmax |delta velocity|={delta['maxVelocityChange']:.3g} m/s\n|delta Cd|={delta['forceXChange']:.3g}",transform=ax.transAxes,fontsize=9)
    fig.suptitle('Native steady Re20: one fixed 32-segment circle and fixed exterior domain',fontsize=15)
    fig.supxlabel('Polygon geometry, finite-domain effects and continuum-reference error remain. This is spatial/iteration sensitivity, not accuracy qualification.',fontsize=10)
    args.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(str(args.output)+'.png',dpi=145);plt.close(fig)
    native.write_json(Path(str(args.output)+'.json'),{**report,'studySha256':native.sha256_file(source),'plotSourceSha256':native.sha256_file(Path(__file__))})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--study',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args())
