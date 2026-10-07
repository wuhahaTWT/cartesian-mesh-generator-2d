#!/usr/bin/env python3
"""Compare native exported states on the same mesh; never assembles a PDE."""
import argparse,csv,json,pathlib

def rows(path):
    with open(path) as stream:return list(csv.DictReader(stream))

def field_difference(left,right):
    a,b=rows(str(left)+'.cells.csv'),rows(str(right)+'.cells.csv')
    if len(a)!=len(b) or any(any(x[k]!=y[k] for k in ('cell','x','y','area')) for x,y in zip(a,b)):
        raise ValueError('field comparison requires identical exported mesh')
    volume=sum(float(x['area']) for x in a)
    return {k:{'volumeL1':sum(float(x['area'])*abs(float(x[k])-float(y[k])) for x,y in zip(a,b))/volume,
               'maximum':max(abs(float(x[k])-float(y[k])) for x,y in zip(a,b))} for k in ('rho','u','v','p','T')}

def compare(candidate,strict,half):
    a,b=rows(str(candidate)+'.history.csv'),rows(str(strict)+'.history.csv')
    if not a or not b or a[-1]['time']!=b[-1]['time']:raise ValueError('terminal clocks differ')
    result={'sameAcceptedClocks':[x['time'] for x in a]==[x['time'] for x in b],
            'candidateMinusStrict':field_difference(candidate,strict)}
    if half:
        h=rows(str(half)+'.history.csv')
        if not h or h[-1]['time']!=b[-1]['time']:raise ValueError('half reference terminal clock differs')
        result.update(strictMinusHalf=field_difference(strict,half),candidateMinusHalf=field_difference(candidate,half))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('candidate',type=pathlib.Path);p.add_argument('strict',type=pathlib.Path);p.add_argument('--half',type=pathlib.Path);a=p.parse_args()
    print(json.dumps(compare(a.candidate,a.strict,a.half),indent=2))
