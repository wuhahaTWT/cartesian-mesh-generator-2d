#!/usr/bin/env python3
"""Summarize native output, not an independent equation/topology audit."""
import argparse, collections, csv, hashlib, json
from pathlib import Path

def rows(path):
    with path.open() as f:return list(csv.DictReader(f))

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main(a):
    root=Path(a.input);cases={};series={}
    for directory in sorted(root.iterdir()):
        path=directory/'result.metrics.json'
        if not path.exists():continue
        meta=json.loads(path.read_text());assert meta['run']['code']==0 and meta['outputIntegrity'],path
        for name,sha in meta['outputSha256'].items():assert digest(directory/name)==sha,(directory,name)
        heat=rows(directory/'result.heat-history.csv');attempts=rows(directory/'result.attempt-history.csv')
        cells=rows(directory/'result.cells.csv');boundary=rows(directory/'result.boundary-heat-history.csv')
        area=sum(float(r['area']) for r in cells);s=meta['summary']
        gain=meta['heatContent']-300*area
        flux_gain=sum(float(r['dt'])*(float(r['sourceIntegral'])-float(r['boundaryFlux'])) for r in heat)
        defect=sum(float(r['dt'])*float(r['globalBalance']) for r in heat)
        wall=[r for r in boundary if r['name'] in ['wall','lid','bottom']]
        end_wall=sum(float(r['diffusiveFlux']) for r in wall if float(r['time'])==s['time'])
        cases[directory.name]={
            'case':meta['case'],'cells':meta['cells'],'endTime':s['time'],'acceptedSteps':s['completedSteps'],
            'rejectedAttempts':s['rejectedAttempts'],'seconds':meta['run']['seconds'],
            'flowSpeed':s['flowSpeed'],'diffusivity':s['diffusivity'],'maxDt':s['maximumTimeStep'],
            'control':s['timeStepControl'],'temperatureScale':s['temperatureScale'],'relativeTimeTolerance':s['timeRelativeTolerance'],
            'heatGain':gain,'integratedSourceMinusOutwardFlux':flux_gain,'integratedNativeBudgetDefect':defect,
            'budgetAccumulationRounding':gain-flux_gain-defect,'endOutwardWallDiffusiveFlux':end_wall,
            'endTemperatureRateIntegral':float(heat[-1]['temporalIntegral']),
            'temperatureRangeOverHistory':[min(float(r['minValue']) for r in heat),max(float(r['maxValue']) for r in heat)],
            'endTemperatureRange':[s['minValue'],s['maxValue']],
            'outletTemperatureRise':meta.get('outletTemperature',300)-300 if 'outletTemperature' in meta else None,
            'rejectionReasons':dict(collections.Counter(r['reason'] for r in attempts if r['reason']!='accepted')),
            'solverQualityPassed':s['solverQualityPassed'],'nativeTopologyRevalidated':s['nativeTopologyRevalidated'],
            'outputIntegrity':True,'binarySha256':meta['provenance']['binarySha256'],
            'checkpointSha256':digest(directory/'result.thermal.checkpoint'),
            'rawDirectory':str(directory)}
        series[directory.name]=(heat,area)
    def field_difference(left,right):
        x=rows(root/left/'result.cells.csv');y=rows(root/right/'result.cells.csv');assert len(x)==len(y)
        return max(abs(float(p['value'])-float(q['value'])) for p,q in zip(x,y))
    e1=field_difference('cavity-dt-1','cavity-dt-05');e2=field_difference('cavity-dt-05','cavity-dt-025')
    comparison={}
    for kind in ['channel','cavity','cylinder']:
        coarse=cases[kind+'-coarse'];fine=cases[kind+'-fine']
        comparison[kind]={'heatGainDifference':fine['heatGain']-coarse['heatGain'],
            'heatGainRelativeDifferenceToFine':abs(fine['heatGain']-coarse['heatGain'])/abs(fine['heatGain']),
            'wallFluxDifference':fine['endOutwardWallDiffusiveFlux']-coarse['endOutwardWallDiffusiveFlux']}
        if coarse['outletTemperatureRise'] is not None:comparison[kind]['outletRiseDifference']=fine['outletTemperatureRise']-coarse['outletTemperatureRise']
    report={'scope':'native 2D incompressible constant-property laminar flow / one-way passive temperature',
        'units':{'temperature':'K','time':'s','heatContent':'K m^2 per unit depth (multiply by rho*cp for energy)',
                 'heatFlux':'K m^2/s per unit depth; outward positive'},
        'cases':cases,'gridComparison':comparison,
        'timeSensitivity':{'case':'cavity','maximumDt':[.1,.05,.025],'eventClipped':True,
            'adjacentEndpointMaxTemperatureDifferences':[e1,e2],'differenceRatio':e1/e2,
            'interpretation':'Consistent with first-order BE; not an adaptive-global-error guarantee.',
            'adaptiveRtol':[.01,.005,.0025],
            'adaptiveAdjacentEndpointMaxDifferences':[field_difference('cavity-coarse','cavity-tight'),field_difference('cavity-tight','cavity-tighter')],
            'adaptiveInterpretation':'Endpoint differences need not decrease monotonically: local control changes the step sequence, and maximum dt / absolute tolerance also limit refinement.',
            'finerReferenceMaximumDt':.0125,
            'adaptiveEndpointDifferencesToFinerReference':[field_difference(k,'cavity-dt-0125') for k in ['cavity-coarse','cavity-tight','cavity-tighter']],
            'finestTwoFixedEndpointDifference':field_difference('cavity-dt-025','cavity-dt-0125')},
        'limitations':['Two spatial levels diagnose sensitivity, not grid independence or general curved-wall heat-transfer accuracy.',
            'Sharp cylinder heating produces a transient spatial undershoot below 300 K; explicitly retained, no clipping. Implicit nonorthogonal corrected diffusion is not a discrete maximum-principle guarantee.',
            'Local BE step-doubling tolerance is not a bound on accumulated, spatial, or unresolved startup error.',
            'Observed live output replacement is detected and isolated; the external process responsible is not established.']}
    output=Path(a.output);output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(report,indent=2)+'\n')
    if a.plot:
        import matplotlib;matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(2,3,figsize=(13,6.4),layout='constrained')
        for col,kind in enumerate(['channel','cavity','cylinder']):
            for level,style in [('coarse','--'),('fine','-')]:
                h,area=series[kind+'-'+level];t=[float(r['time']) for r in h]
                axes[0,col].plot(t,[float(r['heatContent'])/area-300 for r in h],style,label=f"{level}: {cases[kind+'-'+level]['cells']} cells")
                axes[1,col].plot(t,[float(r['minValue'])-300 for r in h],style,color='#2166ac',alpha=.8,label=f'{level} minimum')
                axes[1,col].plot(t,[float(r['maxValue'])-300 for r in h],style,color='#b2182b',alpha=.8,label=f'{level} maximum')
            axes[0,col].set_title(kind.capitalize());axes[0,col].legend(fontsize=8)
            axes[0,col].set_ylabel('Volume mean temperature rise (K)');axes[1,col].set_ylabel('Temperature range relative to 300 K')
            for ax in axes[:,col]:ax.set_xlabel('Physical time (s)');ax.grid(alpha=.2)
        axes[1,0].legend(fontsize=7,loc='center right')
        axes[1,2].text(.03,.12,'Early undershoot: -0.029 / -0.024 K\nRetained; no temperature clipping',transform=axes[1,2].transAxes,fontsize=8)
        fig.suptitle('Native accepted states: heating, cooling and time events')
        fig.savefig(a.plot,dpi=160);plt.close(fig)
    print(json.dumps({'cases':len(cases),'timeSensitivity':report['timeSensitivity'],'gridComparison':comparison},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True);p.add_argument('--plot');main(p.parse_args())
