#!/usr/bin/env python3
"""Summarize exported native fields/fluxes only; no independent PDE operator."""
import argparse,pathlib,csv,json,numpy as np
p=argparse.ArgumentParser();p.add_argument('--root',type=pathlib.Path,required=True);a=p.parse_args();root=a.root
names=['corrected','hybrid-heat','hybrid']
def field(path):
 rows=list(csv.DictReader(path.open()));return {k:np.array([float(r[k]) for r in rows]) for k in rows[0]}
def delta(x,y,ids=None):
 if ids is None:ids=np.arange(len(x['area']))
 w=x['area'][ids];return {k:{'areaL1':float(np.sum(w*np.abs(x[k][ids]-y[k][ids]))/w.sum()),'maximum':float(np.max(np.abs(x[k][ids]-y[k][ids])))} for k in ['rho','u','v','p','T']}
sets={'channel-start':root/'shared-channel'/'start','channel-start-half':root/'shared-channel-half'/'start','channel-full':root/'shared-channel-full'/'comparison','channel-full-half':root/'shared-channel-full-half'/'comparison','cylinder-start':root/'shared-cylinder'/'comparison','cylinder-start-half':root/'shared-cylinder-half'/'comparison','cylinder-full':root/'shared-cylinder-full'/'comparison','cylinder-full-half':root/'shared-cylinder-full-half'/'comparison'}
summary={};fields={}
for label,prefix in sets.items():
 if not prefix.with_suffix('.json').exists():continue
 overall=json.loads(prefix.with_suffix('.json').read_text());report={'overall':overall,'schemes':{}};fields[label]={}
 for name in names:
  stem=pathlib.Path(str(prefix)+'.'+name);x=field(pathlib.Path(str(stem)+'.cells.csv'));fields[label][name]=x;j=json.loads(pathlib.Path(str(stem)+'.json').read_text());history=list(csv.DictReader(pathlib.Path(str(stem)+'.history.csv').open()));bh=list(csv.DictReader(pathlib.Path(str(stem)+'.boundary-history.csv').open()));dt={int(r['step']):float(r['dt']) for r in history};groups={}
  for r in bh:
   g=groups.setdefault(r['name'],{'timeIntegral':{k:0. for k in ['mass','mx','my','energy','heat','viscousWork']},'lastFlux':{k:0. for k in ['mass','mx','my','energy','heat','viscousWork']}})
   for k in g['timeIntegral']:g['timeIntegral'][k]+=dt[int(r['step'])]*float(r[k]);g['lastFlux'][k]+=float(r[k]) if int(r['step'])==int(history[-1]['step']) else 0.
  j['boundaryGroups']=groups;j['rangesLastAccepted']={k:[float(x[k].min()),float(x[k].max())] for k in ['rho','p','T']};j['rangesAllAccepted']={k:[min(float(r[k+'Min']) for r in history),max(float(r[k+'Max']) for r in history)] for k in ['rho','p','t']} if history else {};j['acceptedTimes']=[float(r['time']) for r in history];j['maximumAbsoluteGlobalBalance']={k:max(abs(float(r[k])) for r in history) for k in ['balanceMass','balanceMx','balanceMy','balanceEnergy']} if history else {};j['maximumCellBalance']=max((float(r['maxCellBalance']) for r in history),default=0.);report['schemes'][name]=j
  if pathlib.Path(str(stem)+'.faces.csv').exists():
   faces=list(csv.DictReader(pathlib.Path(str(stem)+'.faces.csv').open()));
   # Boundary grouping is authoritative for force/heat; near-wall IDs use
   # prescribed no-slip face IDs from the original boundary file, not this heuristic.
   bpath=root/'inputs'/('channel.boundaries' if label.startswith('channel') else 'cylinder.boundaries');wallIds={int(r.split()[0]) for r in bpath.read_text().splitlines() if ' no-slip-wall ' in r};ids=np.array(sorted({int(r['owner']) for r in faces if int(r['face']) in wallIds}),dtype=int);w=x['area'][ids];j['nearWall']={k:{'areaMean':float(np.sum(w*x[k][ids])/w.sum()),'minimum':float(x[k][ids].min()),'maximum':float(x[k][ids].max())} for k in ['T','p']};j['nearWallCells']=ids.tolist()
 assert report['schemes'][names[0]]['acceptedTimes']==report['schemes'][names[1]]['acceptedTimes']==report['schemes'][names[2]]['acceptedTimes'];report['actualClockMatches']=True
 for name in names[1:]:report['schemes'][name]['differenceFromCorrected']=delta(fields[label][name],fields[label]['corrected']);report['schemes'][name]['nearWallDifferenceFromCorrected']=delta(fields[label][name],fields[label]['corrected'],np.array(report['schemes'][name].get('nearWallCells',np.arange(len(fields[label][name]['area'])))))
 summary[label]=report
for label in ['channel-start','channel-full','cylinder-start','cylinder-full']:
 if label in fields and label+'-half' in fields and summary[label]['overall']['status']=='complete' and summary[label+'-half']['overall']['status']=='complete' and summary[label]['overall']['physicalTime']==summary[label+'-half']['overall']['physicalTime']:
  summary[label]['temporalHalfStepDifference']={n:delta(fields[label][n],fields[label+'-half'][n]) for n in names}
(root/'comparison-summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps({label:{n:{'ranges':r['rangesLastAccepted'],'seconds':r['exclusiveSecondsWithExport'],'discarded':r['commonTransactionDiscarded']} for n,r in v['schemes'].items()} for label,v in summary.items()},indent=2))
