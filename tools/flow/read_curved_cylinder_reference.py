#!/usr/bin/env python3
"""Read actual native geometry/flux exports; no independent PDE discretization."""
import pathlib,json,csv,numpy as np,math,gzip
root=pathlib.Path('outputs/corrected-curved-reference');old=pathlib.Path('outputs/curved-cylinder')
rho=1.176624281484062;nu=1.846e-5/rho;alpha=.025759/(rho*1.4*287.05/.4)
def rows(p):return list(csv.DictReader(p.open()))
def wstats(x,w):
 o=np.argsort(x);c=np.cumsum(w[o])/sum(w);return {'minimum':float(min(x)),'p50':float(x[o[np.searchsorted(c,.5)]]),'p95':float(x[o[np.searchsorted(c,.95)]]),'maximum':float(max(x))}
def read(prefix):
 j=json.loads(prefix.with_suffix('.json').read_text());cs=rows(prefix.with_suffix('.cells.csv'));fs=rows(prefix.with_suffix('.faces.csv'));walls=[f for f in fs if f.get('kind')=='no-slip-wall'];length=np.array([math.hypot(float(f['sx']),float(f['sy'])) for f in walls]);d=[];normal=[];skew=[];temp=[];pressure=[];heat=[];shear=[]
 geo=json.loads((prefix.parent/'mesh.resolution.json').read_text()) if (prefix.parent/'mesh.resolution.json').exists() else json.loads((old/('N192-W48' if prefix.name.startswith('mid') else 'N192-W96')/'mesh.resolution.json').read_text())
 for f,l in zip(walls,length):
  c=cs[int(f['owner'])];n=np.array([float(f['sx']),float(f['sy'])])/l;v=np.array([float(f['x'])-float(c['x']),float(f['y'])-float(c['y'])]);dn=float(np.dot(n,v));d.append(dn);skew.append(float(np.linalg.norm(v-n*dn)/dn));temp.append(float(c['temperature']));pressure.append(float(c['p']));heat.append(float(f['heatFlux']));t=np.array([-n[1],n[0]]);shear.append(float(np.dot(t,[float(f['viscousMomentumX']),float(f['viscousMomentumY'])])/l))
 hp=prefix.with_suffix('.history.csv');hz=prefix.with_suffix('.history.csv.gz');history=list(csv.DictReader(hp.open())) if hp.exists() else list(csv.DictReader(gzip.open(hz,'rt')));balances={k:max(abs(float(h[k])) for h in history) for k in ['cellBalanceError','balanceMass','balanceMomentumX','balanceMomentumY','balanceEnergy']}
 time=j['time'];heat=np.array(heat);shear=np.array(shear);deltaT=math.sqrt(alpha*time);deltaV=math.sqrt(nu*time)
 result={'time':time,'status':j['status'],'failure':j['failure'],'targetReached':j['targetReached'],'cells':j['cells'],'requestedMaximumStep':j['maximumStep'],'initialTime':j['initialTime'],'acceptedSteps':j['acceptedSteps'],'rejected':j['rejectedCandidates'],'solverSeconds':j['elapsedSeconds'],'spatialCalls':j['spatialEvaluations'],'maximumAbsoluteNativeBalances':balances, 'rhoMin':min(float(c['rho']) for c in cs),'pMin':min(float(c['p']) for c in cs),'temperatureRange':[min(float(c['temperature']) for c in cs),max(float(c['temperature']) for c in cs)],'nearWallTemperature':wstats(np.array(temp),length),'centroidNormalDistanceMeters':wstats(np.array(d),length),'normalNonorthogonalityTangentOverNormal':wstats(np.array(skew),length),'diffusionScalesMeters':{'sqrtNuT':deltaV,'sqrtAlphaT':deltaT},'resolutionNormalExtentMeters':{k:v*geo['reference_length'] for k,v in geo['actual']['wall_owner_normal_extent_over_reference'].items()},'wallHeat':{'netWperMeter':float(sum(heat)),'absoluteSumWperMeter':float(sum(abs(heat))),'positiveWperMeter':float(sum(heat[heat>0])),'negativeWperMeter':float(sum(heat[heat<0]))},'wallShearPa':wstats(shear,length),'wallForceNperMeter':{k:sum(float(f[k]) for f in walls) for k in ['momentumX','momentumY','viscousMomentumX','viscousMomentumY']},'wallLengthMeters':float(sum(length))}
 # Net heat cancellation and angular nonuniformity are export diagnostics only.
 bins=[]
 for lo in np.linspace(-math.pi,math.pi,33)[:-1]:
  ids=[i for i,f in enumerate(walls) if lo<=math.atan2(float(f['y']),float(f['x']))<lo+math.pi/16]
  bins.append({'theta':lo+math.pi/32,'length':float(sum(length[ids])),'heat':float(sum(heat[ids])),'meanTemperature':float(np.dot(length[ids],np.array(temp)[ids])/sum(length[ids])) if ids else None,'shearIntegral':float(np.dot(length[ids],shear[ids]))})
 result['angular32']=bins
 return result
out={'scope':'Geometry and native export readback, not PDE reconstruction or independent physical reference','radiusMeters':5e-5,'nu':nu,'alpha':alpha,'cases':{}}
for name in ['N192-W16/linear','N192-W48/linear-resume','N192-W96/linear','N192-W96/linear-half','N48-W48/linear']:
 out['cases'][name]=read(old/name)
for name in ['mid-resumed','mid-1ns','fine','fine-half']:
 p=root/name
 if p.with_suffix('.json').exists():out['cases']['new/'+name]=read(p)
out['chordSagittaMeters']={str(n):5e-5*(1-math.cos(math.pi/n)) for n in [48,192]}
for label,a,b in [('old-space','N192-W48/linear-resume','N192-W96/linear'),('old-time','N192-W96/linear','N192-W96/linear-half'),('new-space','new/mid-1ns','new/fine'),('new-time','new/fine','new/fine-half')]:
 if a in out['cases'] and b in out['cases']:
  x,y=out['cases'][a],out['cases'][b];out[label]={'sameEndTime':x['time']==y['time'],'bothReached':x['targetReached'] and y['targetReached'],'netHeatDifference':x['wallHeat']['netWperMeter']-y['wallHeat']['netWperMeter'],'relativeNetHeatDifference':abs(x['wallHeat']['netWperMeter']-y['wallHeat']['netWperMeter'])/abs(y['wallHeat']['netWperMeter']),'angularHeatAbsoluteDifference':sum(abs(i['heat']-j['heat']) for i,j in zip(x['angular32'],y['angular32'])),'forceXDifference':x['wallForceNperMeter']['momentumX']-y['wallForceNperMeter']['momentumX'],'viscousForceXDifference':x['wallForceNperMeter']['viscousMomentumX']-y['wallForceNperMeter']['viscousMomentumX']}
(root/'readback.json').write_text(json.dumps(out,indent=2))
print(json.dumps({n:{k:v[k] for k in ['time','cells','wallHeat','wallForceNperMeter','centroidNormalDistanceMeters','normalNonorthogonalityTangentOverNormal']} for n,v in out['cases'].items()},indent=2))
