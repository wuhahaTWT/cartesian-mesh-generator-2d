#!/usr/bin/env python3
"""Measure audited flat-plate fields; reference differences are NOT acceptance gates.

Station profiles are piecewise-linear in x between cell centres, with no
extrapolation. This diagnostic deliberately supports complete axis-aligned
rectilinear meshes only. It does not sample arbitrary Cut-cell meshes.
"""
import argparse
import bisect
import hashlib
import json
import math
import re
from pathlib import Path

import verify_sst_rans as sst

STATIONS = (.97008, 1.90334)
REFERENCE_BASE = 'https://tmbwg.github.io/turbmodels/'
REFERENCE_FILES = {
    'cf_plate_sstv.dat': 'FlatPlate/SST/',
    'flatplate_u_sstv.dat': 'FlatPlate/SST/',
    'cf_incomp_results_sstv.dat': 'FlatPlate_validation/',
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def finite(values):
    return all(math.isfinite(v) for v in values)


def interpolate(x, y, target):
    require(len(x) == len(y) and len(x) >= 2 and finite([*x, *y, target]), 'invalid interpolation data')
    require(all(b > a for a, b in zip(x, x[1:])), 'interpolation coordinates must increase')
    require(x[0] <= target <= x[-1], 'profile extrapolation is not supported')
    j = bisect.bisect_left(x, target)
    if j < len(x) and x[j] == target:
        return y[j]
    f = (target - x[j-1]) / (x[j] - x[j-1])
    return y[j-1] + f * (y[j] - y[j-1])


def first_crossing(y, u, target=.99):
    require(len(y) == len(u) and len(y) >= 2 and target > 0 and
            finite([*y, *u, target]) and y[0] == 0 and u[0] == 0,
            'wall-anchored finite profile required')
    require(all(b > a for a, b in zip(y, y[1:])), 'profile heights must increase')
    for j in range(1, len(y)):
        if u[j-1] < target <= u[j]:
            return y[j-1] + (y[j] - y[j-1]) * (target-u[j-1]) / (u[j]-u[j-1])
    return None


def tensor_field(mesh, measured, rows):
    xs = sorted({v[0] for v in mesh.vertices})
    ys = sorted({v[1] for v in mesh.vertices})
    nx, ny = len(xs)-1, len(ys)-1
    require(nx*ny == len(mesh.cells) == len(rows), 'complete tensor mesh required')
    xi, yi = {x:i for i,x in enumerate(xs)}, {y:i for i,y in enumerate(ys)}
    grid = {}
    xc, yc = {}, {}
    for cell, centre, row in zip(mesh.cells, measured.centroids, rows):
        verts = [mesh.vertices[i] for i in cell.vertices]
        left, right = min(v[0] for v in verts), max(v[0] for v in verts)
        bottom, top = min(v[1] for v in verts), max(v[1] for v in verts)
        require(len(verts) == 4 and set(verts) == {(left,bottom),(right,bottom),(right,top),(left,top)},
                'axis-aligned quadrilateral required')
        i, j = xi[left], yi[bottom]
        require(xi[right] == i+1 and yi[top] == j+1 and (i,j) not in grid, 'non-tensor or duplicate cell')
        require(int(row['cell']) == cell.id, 'cell ordering differs')
        require(xc.get(i, centre[0]) == centre[0] and yc.get(j, centre[1]) == centre[1], 'misaligned tensor centres')
        xc[i], yc[j] = centre
        grid[i,j] = {key:sst.num(row[key],key) for key in ('u','v','p','k','omega','nuT')}
    require(len(grid) == nx*ny, 'incomplete tensor field')
    return [xc[i] for i in range(nx)], [yc[j] for j in range(ny)], grid


def station_profile(xc, yc, grid, wall_samples, station, nu, speed, wall_y=0):
    require(nu > 0 and speed > 0 and finite([nu,speed,wall_y]), 'invalid profile scales')
    require(yc[0] > wall_y, 'cell centres must be above wall')
    samples = sorted(wall_samples, key=lambda r:r['x'])
    tau = interpolate([r['x'] for r in samples], [r['kinematicShear'] for r in samples], station)
    require(tau > 0, 'wall scaling requires positive shear; separated profiles unsupported')
    utau = math.sqrt(tau)
    values = []
    for j, y in enumerate(yc):
        fields = {key:interpolate(xc, [grid[i,j][key] for i in range(len(xc))], station)
                  for key in ('u','v','p','k','omega','nuT')}
        values.append({'y':y-wall_y, **fields, 'uOverU':fields['u']/speed,
                       'uPlus':fields['u']/utau, 'yPlus':(y-wall_y)*utau/nu})
    return {'x':station, 'Cf':2*tau/speed**2, 'uTau':utau,
            'delta99FirstCrossing':first_crossing([0]+[v['y'] for v in values], [0]+[v['uOverU'] for v in values]),
            'samples':values}


def file_record(path):
    path = Path(path)
    return {'path':str(path), 'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}


def analyze_case(name, mesh_path, prefix, stations=STATIONS):
    # Recompute every original-equation gate from actual input files. A saved
    # report's "valid" flag is never sufficient to unlock physical plotting.
    audit = sst.audit(mesh_path, prefix)
    mesh, measured, meta, rows, _, _ = sst.read_artifacts(mesh_path, prefix)
    require(meta['case'] == 'flatplate', 'flat plate case required')
    require(all(x > meta['flatPlateLeadingEdge'] for x in stations), 'station is not on the plate')
    xc, yc, grid = tensor_field(mesh, measured, rows)
    profiles = [station_profile(xc,yc,grid,audit['plateWallSamples'],x,meta['nu'],meta['speed'],measured.bounds[1])
                for x in stations]
    wall_length = sum(r['length'] for r in audit['plateWallSamples'])
    drag = sum(r['kinematicShear']*r['length'] for r in audit['plateWallSamples'])
    files = [file_record(mesh_path)] + [file_record(str(prefix)+suffix)
             for suffix in ('.json','.cells.csv','.faces.csv','.history.csv')]
    return {'name':name, 'equationAuditValid':True, 'physicalAccuracyQualified':False,
            'cells':len(mesh.cells), 'nx':len(xc), 'ny':len(yc), 'bounds':list(measured.bounds),
            'physicalInputs':audit['physicalInputs'], 'model':meta['model'],
            'leadingEdge':meta['flatPlateLeadingEdge'], 'topBoundary':meta['flatPlateTop'],
            'convection':meta['convection'], 'wallLength':wall_length,
            'CdSkinFriction':2*drag/(meta['speed']**2*wall_length),
            'wallSamples':audit['plateWallSamples'], 'profiles':profiles, 'files':files,
            'audit':audit}


def read_zones(path, columns):
    """Read the small ASCII POINT-zone tables published for this TMR case."""
    zones, current = {}, None
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith('#') or line.lower().startswith('variables'):
            continue
        if line.lower().startswith('zone'):
            match = re.search(r't\s*=\s*"([^"]+)"', line, re.I)
            require(match is not None and match[1] not in zones, 'invalid or duplicate reference zone')
            current = match[1]
            zones[current] = []
        else:
            require(current is not None, 'reference values before zone')
            values = [float(v) for v in line.split()]
            require(len(values) == columns and finite(values), 'invalid reference row')
            zones[current].append(values)
    require(zones and all(len(v) >= 2 for v in zones.values()), 'empty reference zone')
    return zones


def load_references(folder):
    zones = {name:read_zones(folder/name, 4 if 'incomp' in name else 2) for name in REFERENCE_FILES}
    return {'sourcePage':REFERENCE_BASE+'flatplate_sst.html',
            'scope':'Context only: SST/SST-V/SST-Vm, not a matched SST-2003m verification target',
            'compressibleMach':.2, 'lengthScale':1., 'reynoldsPerUnitLength':5e6,
            'velocityTableScale':'Interpreted as u/U_infinity and y/L, L=1, from the nondimensional case and far-field values near 1; the table itself labels only u,y',
            'files':[{**file_record(folder/name), 'url':REFERENCE_BASE+subdir+name} for name,subdir in REFERENCE_FILES.items()],
            'tables':zones}


def reference_comparison(case, reference):
    require(case['leadingEdge'] == 0 and case['bounds'][1] == 0 and case['wallLength'] == 2,
            'reference coordinate convention differs')
    p = case['physicalInputs']
    require(abs(p['speed']/p['nu']/5e6-1) < 1e-12 and
            abs(p['inletK']/p['speed']**2/2.25e-7-1) < 1e-12 and
            abs(p['inletOmega']/p['speed']/125-1) < 1e-12, 'reference physical scale differs')
    comparisons = []
    for profile in case['profiles']:
        station = profile['x']
        cf = {}
        for name, pairs in reference['tables']['cf_plate_sstv.dat'].items():
            expected = interpolate([v[0] for v in pairs], [v[1] for v in pairs], station)
            cf[name] = {'reference':expected, 'signedRelativeDifference':profile['Cf']/expected-1}
        profile_tables = reference['tables']['flatplate_u_sstv.dat']
        key = next((key for key in profile_tables if float(key.split('=')[1]) == station), None)
        require(key is not None, 'no reference station at requested coordinate')
        pairs = profile_tables[key]
        samples = [{ 'y':s['y'], 'nativeUOverU':s['uOverU'],
                    'referenceUOverU':interpolate([v[1] for v in pairs],[v[0] for v in pairs],s['y'])}
                   for s in profile['samples'] if 0 < s['y'] <= .05]
        require(samples, 'no near-wall reference overlap')
        comparisons.append({'x':station,'cf':cf,'profileSamples':samples,
                            'maxAbsoluteVelocityDifferenceOverU':max(abs(s['nativeUOverU']-s['referenceUOverU']) for s in samples)})
    return {'qualified':False, 'reason':'Different SST model variants, compressibility, numerical grids and incompletely matched boundary setup; descriptive differences only',
            'stations':comparisons}


def compare_cases(a, b):
    for key in ('physicalInputs','model','leadingEdge','topBoundary','convection'):
        require(a[key] == b[key], 'case physics differs: '+key)
    require([p['x'] for p in a['profiles']] == [p['x'] for p in b['profiles']], 'station sets differ')
    result = []
    for pa,pb in zip(a['profiles'],b['profiles']):
        sa,sb = pa['samples'],pb['samples']
        lower = max(sa[0]['y'],sb[0]['y'])
        upper = min(.05,sa[-1]['y'],sb[-1]['y'])
        require(0 < lower < upper, 'no shared near-wall sample interval')
        ys = [lower*(upper/lower)**(i/99) for i in range(100)]
        ys[0],ys[-1] = lower,upper
        differences = [interpolate([v['y'] for v in sb],[v['uOverU'] for v in sb],y)-
                       interpolate([v['y'] for v in sa],[v['uOverU'] for v in sa],y) for y in ys]
        result.append({'x':pa['x'],'CfA':pa['Cf'],'CfB':pb['Cf'],
                       'signedRelativeCfChange':pb['Cf']/pa['Cf']-1,
                       'profileCommonY':ys, 'profileDifferenceOverU':differences,
                       'maxAbsoluteVelocityChangeOverU':max(map(abs,differences))})
    return {'from':a['name'],'to':b['name'],'stations':result,
            'signedRelativeCdChange':b['CdSkinFriction']/a['CdSkinFriction']-1,
            'interpretation':'Sensitivity observation, not a convergence order or error estimate'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--reference-dir',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    names = [c['name'] for c in manifest['cases']]
    require(names and len(set(names)) == len(names), 'case names must be unique')
    cases = [analyze_case(c['name'],Path(c['mesh']),Path(c['prefix'])) for c in manifest['cases']]
    by_name = {c['name']:c for c in cases}
    result = {'scope':'Audited rectilinear SST flat-plate physical diagnostics',
              'physicalAccuracyQualified':False, 'manifest':file_record(args.manifest), 'cases':cases,
              'comparisons':[compare_cases(by_name[a],by_name[b]) for a,b in manifest.get('comparisons',[])],
              'limits':['No extrapolation; cell-centre profiles interpolated linearly in x',
                        'delta99 is the first crossing of 0.99 U_infinity, not a measured edge velocity',
                        'Reference differences do not set or relax any solver acceptance gate',
                        'No grid convergence order inferred from these non-similar meshes']}
    if args.reference_dir:
        result['reference'] = load_references(args.reference_dir)
        for c in cases:
            c['referenceComparison'] = reference_comparison(c,result['reference'])
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'cases':[{k:c[k] for k in ('name','cells','CdSkinFriction')} for c in cases],
                      'physicalAccuracyQualified':False}))


if __name__ == '__main__':
    main()
