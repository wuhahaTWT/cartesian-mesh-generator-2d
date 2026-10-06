"""Real Cut-cell native regression. Input preparation and output reads only."""
import json,math,pathlib,subprocess,sys,tempfile,csv
repo,build=map(pathlib.Path,sys.argv[1:3])
with tempfile.TemporaryDirectory(prefix='cartmesh-hmm-') as d:
 r=pathlib.Path(d);xy=r/'ring.xy'
 xy.write_text('\n\n'.join('\n'.join(f'{radius*math.cos(2*math.pi*(i+.17)/128):.17g} {radius*math.sin(2*math.pi*(i+.17)/128):.17g}' for i in range(128)) for radius in [1.,.5])+'\n')
 mesh=r/'mesh'
 subprocess.run([str(build/'cartmesh2d_cli'),str(xy),str(mesh),'4','.125','.1','interior',str(r/'foam'),'3','0'],check=True,timeout=30)
 prefix=r/'hmm';subprocess.run([str(build/'cartmesh2d_hybrid_heat_benchmark'),str(mesh)+'.solver.cm2d',str(prefix),'trace','.05','.005'],check=True,timeout=30)
 j=json.loads(prefix.with_suffix('.json').read_text())
 assert j['cells']==264 and j['steps']==10 and j['time']==.05
 assert j['finalPerturbationEnergy']<j['initialPerturbationEnergy'] and j['maximumPerturbationEnergyIncrease']==0
 assert j['minimumTemperatureAllStates']>0
 # Trace residual is in W/m per depth, same dimensional budget as reused PCG.
 assert j['steadyFaceFluxJump']<4*(1e-13+1e-12*j['traceRhsNorm'])
 # This mesh/k/time has a resolved spatial error: solver tolerance must not
 # contaminate that measurement. No universal spatial-accuracy threshold.
 assert j['toleranceSensitivityK']<.01*.2*j['temperatureRelativeL1']
 rows=list(csv.DictReader(prefix.with_suffix('.history.csv').open()))
 assert len(rows)==10
 # Triangle bound from reported native cell residuals and face jumps, W/m.
 for x in rows:
  bound=j['cells']*float(x['cellBalance'])+j['traceUnknowns']*float(x['faceFluxJump'])+32768*sys.float_info.epsilon*2*.6707893009904248
  assert abs(float(x['globalHeatBalance']))<=bound
 probe=r/'probe';subprocess.run([str(build/'cartmesh2d_hybrid_heat_benchmark'),str(mesh)+'.solver.cm2d',str(probe),'uniform','1e-5','1e-6','impulse=88'],check=True,timeout=30)
 q=json.loads(probe.with_suffix('.json').read_text());assert q['steps']==10 and q['time']==1e-5
 assert q['maximumPrincipleUndershootK']>0 and q['minimumTemperatureAllStates']>0
 # Known non-monotone operator: retain the measured limitation; no clipping.
print('native actual 264-cell HMM geometry/steady/decay/clock regression passed; Euler coupling excluded')
