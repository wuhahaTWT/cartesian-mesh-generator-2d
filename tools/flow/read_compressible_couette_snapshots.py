#!/usr/bin/env python3
"""Read native face/residual exports; no reconstruction or PDE discretization."""
import argparse,csv,json,math,pathlib
p=argparse.ArgumentParser();p.add_argument('prefix');a=p.parse_args();prefix=pathlib.Path(a.prefix)
def rows(suffix):
    with pathlib.Path(str(prefix)+suffix).open() as f:return list(csv.DictReader(f))
reference=rows('.reference-faces.csv'); cells=rows('.cells.csv');area=sum(float(c['area']) for c in cells);report={}
for label in ['point','mean','input','final']:
    native=rows('.snapshot-'+label+'.csv');assert len(native)==len(reference)
    terms={name:[[0.,0.] for _ in cells] for name in ['inviscid','viscous','referencePressure','referenceAdvection','referenceViscous','numericalPressure','numericalAdvectiveRemainder','pressureError','advectiveRemainderError','wallInvError','internalInvError','wallViscError','internalViscError']}
    wallTorque=[0.,0.];pressureTorque=0.;wallMass=0.
    for x,y in zip(native,reference):
        assert x['face']==y['face']
        inv=[float(x[k])-float(x[v]) for k,v in [('mx','viscX'),('my','viscY')]]
        visc=[float(x[k]) for k in ['viscX','viscY']];pre=[float(y[k]) for k in ['pressureX','pressureY']];adv=[float(y[k]) for k in ['advX','advY']];vr=[float(y[k]) for k in ['viscX','viscY']]
        numpre=[float(x[k]) for k in ['pressureX','pressureY']]
        numadv=[a-b for a,b in zip(inv,numpre)]
        scope='wall' if x['neighbour']=='-1' else 'internal'
        vectors=dict(inviscid=inv,viscous=visc,referencePressure=pre,referenceAdvection=adv,referenceViscous=vr)
        vectors.update(numericalPressure=numpre,numericalAdvectiveRemainder=numadv,pressureError=[a-b for a,b in zip(numpre,pre)],advectiveRemainderError=[a-b for a,b in zip(numadv,adv)])
        vectors[scope+'InvError']=[a-b-c for a,b,c in zip(inv,pre,adv)];vectors[scope+'ViscError']=[a-b for a,b in zip(visc,vr)]
        for name,vec in vectors.items():
            for c,sign in [(int(x['owner']),1),(int(x['neighbour']),-1)]:
                if c>=0:
                    for k in range(2):terms[name][c][k]+=sign*vec[k]
        if scope=='wall':
            cx,cy=float(x['x']),float(x['y']);side=0 if cx*cx+cy*cy<2.25 else 1
            wallTorque[side]+=cx*visc[1]-cy*visc[0];pressureTorque+=cx*inv[1]-cy*inv[0];wallMass+=abs(float(y['mass']))
    report[label]=dict(volumeL1Divergence={name:sum(math.hypot(*v) for v in values)/area for name,values in terms.items()},viscousWallTorque=wallTorque,wallInviscidTorque=pressureTorque,continuousReferenceAbsoluteWallMassFlux=wallMass)
print(json.dumps(report,indent=2))
