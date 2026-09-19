#!/usr/bin/env python3
"""Audit fixed-mesh thermal-vortex time refinement against continuous solutions.

This is a finite-time convergence diagnosis, not a universal accuracy threshold.
All cases must retain the same mesh, schemes, properties, tolerances and end time.
"""
import argparse
import json
import math
from pathlib import Path
import verify_thermal_scale as scale


def time_series_checks(cases):
    require=scale.scalar.require
    require(len(cases)>=3,'time study requires at least three stepsizes')
    require(all(c.get('valid') is True for c in cases),'invalid physical audit in time series')
    rows=[dict(c['reference']) for c in cases]
    constants=('meshSha256','cells','h','time','nu','diffusivity','speed','velocityRelaxation',
               'flowConvection','scalarConvection','viscousStress','pressurePreconditioner',
               'coupledFlowTolerance','flowTolerance','scalarRelativeTolerance',
               'scalarAbsoluteTolerance','scalarCellTolerance')
    require(all(all(k in r for k in constants) for r in rows),'time-series metadata missing')
    for row in rows:
        require(row.get('coupledFlowToleranceStatus')=='explicit-thermal-flow-json' and
                row.get('flowToleranceStatus')=='explicit-flow-json','unknown flow tolerance')
        require(row.get('carrierCheckpointByteIdentical') is True,'carrier equivalence missing')
        require(isinstance(row['steps'],int) and not isinstance(row['steps'],bool) and row['steps']>0,
                'invalid accepted step count')
        require(math.isfinite(row['dt']) and row['dt']>0,'invalid time step')
        for key in ('cells','h','time','nu','diffusivity','speed','velocityRelaxation',
                    'coupledFlowTolerance','flowTolerance','scalarRelativeTolerance',
                    'scalarAbsoluteTolerance','scalarCellTolerance'):
            require(math.isfinite(row[key]) and row[key]>0,'invalid time-series '+key)
        scale.scalar.same(row['coupledFlowTolerance'],row['flowTolerance'],'carrier tolerance mismatch',absolute=0,relative=1e-14)
        scale.scalar.same(row['dt']*row['steps'],row['time'],'mixed final time',absolute=1e-14,relative=1e-13)
        for key in constants:
            if key=='time':
                scale.scalar.same(row[key],rows[0][key],'mixed time-study '+key,absolute=1e-14,relative=1e-13)
            else:require(row[key]==rows[0][key],'mixed time-study '+key)
        for key in ('scalarContinuousL2','velocityContinuousL2'):
            require(math.isfinite(row[key]) and row[key]>=0,'invalid continuous-reference error')
    decreasing={key:True for key in ('scalarContinuousL2','velocityContinuousL2')}
    for coarse,fine in zip(rows,rows[1:]):
        require(fine['dt']<coarse['dt'] and fine['steps']>coarse['steps'],'time steps must strictly refine')
        for key in decreasing:
            a,b=coarse[key],fine[key]
            decreasing[key] &= b<a
            fine[key+'Order']=math.log(a/b)/math.log(coarse['dt']/fine['dt']) if a>0 and b>0 else None
    return {'valid':all(decreasing.values()),'decreases':decreasing,'series':rows,
            'scope':'Fixed-mesh continuous-time reference errors must decrease. Observed slopes include spatial and iterative error; no universal temporal-order or engineering-accuracy claim.'}


def audit_pairs(pairs,output):
    report={'valid':False,'cases':[],'issues':[],
            'verifierSha256':scale.native.sha256_file(Path(__file__).resolve())}
    try:
        for index,(prefix,flow) in enumerate(pairs):
            case=scale.audit(prefix,flow,output.parent/f'time-case-{index+1}'/'audit.json')
            report['cases'].append(case)
            scale.scalar.require(case['valid'],'independent thermal/flow audit failed')
            info=scale.read(prefix,'.json');control=scale.read(flow,'.json')
            case['reference'].update(meshSha256=scale.native.sha256_file(Path(info['mesh'])),
                flowConvection=control['convection'],scalarConvection=info['convection'],
                viscousStress=control['viscousStress'],pressurePreconditioner=control['pressurePreconditioner'])
        report['temporal']=time_series_checks(report['cases'])
        report['valid']=report['temporal']['valid']
        if not report['valid']:report['issues'].append('continuous-reference error did not decrease')
    except (ValueError,OSError,KeyError,OverflowError,TypeError) as exc:
        report['issues'].append(str(exc))
    scale.native.write_json(output,report)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case',type=Path,nargs=2,action='append',required=True,metavar=('THERMAL_PREFIX','FLOW_PREFIX'))
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if len(a.case)<3:p.error('at least three --case pairs required')
    r=audit_pairs(a.case,a.output)
    print(json.dumps({'valid':r['valid'],'issues':r['issues']}))
    return 0 if r['valid'] else 1


if __name__=='__main__':raise SystemExit(main())
