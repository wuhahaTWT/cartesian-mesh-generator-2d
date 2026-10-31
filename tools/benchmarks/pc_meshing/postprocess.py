"""Derived diagnostics on retained accepted fields; never solves FV equations."""
from pathlib import Path
import csv,json,math,re,sys
from study import ROOT,BASE,read,dump,sha
sys.path.insert(0,str(ROOT/'tools/flow'))
sys.path.insert(0,str(ROOT/'tools/benchmarks'))
import native_mesh as nm
from annulus_cross import errors,collect_external
import annulus_cross

def fluent_fields(folder):
    record=read(folder/'result.json')
    if record.get('status')!='finished' or not (folder/'cells.csv').exists():return
    mesh=nm.read_cm2d(Path(record['mesh']));g=nm.measure(mesh)
    with (folder/'cells.csv').open(encoding='utf-8-sig') as f:
        source=[{k.strip():v.strip() for k,v in row.items()} for row in csv.DictReader(f)]
    if len(source)!=len(mesh.cells):raise ValueError('Fluent volume export count mismatch')
    if sorted(int(r['cellnumber']) for r in source)!=list(range(1,len(mesh.cells)+1)):
        raise ValueError('Fluent export cell IDs differ from input')
    from fluent_topology import case_cell_map
    topology=case_cell_map(mesh,folder/'accepted.cas.h5')
    import h5py
    with h5py.File(folder/'accepted.dat.h5','r') as f:
        linked=f['settings/Case File'][0].decode()
        if Path(linked).name!='accepted.cas.h5':raise ValueError('data/case pairing changed')
        arrays=[]
        for name in ['SV_U','SV_V','SV_P']:
            d=f['results/1/phase-1/cells/'+name+'/1']
            if int(d.attrs['minId'][0])!=1 or int(d.attrs['maxId'][0])!=len(mesh.cells):
                raise ValueError('unexpected Fluent field index range')
            arrays.append(d[:])
    velocity_key=lambda u,v:(format(float(u),'.9E'),format(float(v),'.9E'))
    exported={velocity_key(row['x-velocity'],row['y-velocity']):row for row in source}
    saved_keys=[velocity_key(u,v) for u,v in zip(arrays[0],arrays[1])]
    if len(exported)!=len(source) or len(set(saved_keys))!=len(source) or set(saved_keys)!=set(exported):
        record.pop('circleErrors',None);record.update(finiteWholeVolumeField=False,cellMappingError='accepted data and ASCII velocity fields do not have a unique bijection')
        dump(folder/'result.json',record);return
    rows=[];offset=0.;same_ids=0;pressure_offsets=[]
    for data_id,(u,v,pressure,k) in enumerate(zip(*arrays,saved_keys),1):
        i=topology[data_id];c=g.centroids[i];ascii=exported[k]
        x=float(ascii['x-coordinate']);y=float(ascii['y-coordinate'])
        offset=max(offset,abs(x-c[0]),abs(y-c[1]))
        same_ids+=int(data_id-1==i)
        pressure_offsets.append(float(ascii['pressure'])-float(pressure))
        values=dict(cell=i,x=c[0],y=c[1],area=g.areas[i],u=float(u),v=float(v),p=float(pressure))
        if not all(math.isfinite(value) for value in values.values()):raise ValueError('nonfinite Fluent field')
        rows.append(values)
    rows.sort(key=lambda row:row['cell']);record.pop('cellMappingError',None)
    record.update(circleErrors=errors(rows),finiteWholeVolumeField=True,
        cellMapping=dict(cells=len(rows),exactNodeAndCellIncidenceBijection=True,caseIdsPreserved=same_ids==len(rows),
            pairedCellDataRangeVerified=True,asciiVelocityBijection=True,maxAsciiCentreDifferenceMetres=offset,
            asciiMinusSavedPressureOffsetMin=min(pressure_offsets),asciiMinusSavedPressureOffsetMax=max(pressure_offsets),
            strategy='Paired accepted case/data global IDs, exact vertices and full cell-face incidence; ASCII row IDs are not cell identities; comparisons remove pressure gauge',
            acceptedCaseSha256=sha(folder/'accepted.cas.h5'),acceptedDataSha256=sha(folder/'accepted.dat.h5')))
    dump(folder/'result.json',record)
    with (folder/'comparison.cells.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=['cell','x','y','area','u','v','p']);w.writeheader();w.writerows(rows)

def foam_resources(folder):
    if not (folder/'inner-time.json').exists():return
    text=(folder/'inner-time.json').read_text(encoding='utf-8')
    t={key:float(value) for key,value in re.findall(r'"(\w+)":([0-9]+(?:\.[0-9]+)?)\s*[,}]',text)}
    status=read(folder/'inner-status.json') if (folder/'inner-status.json').exists() else {}
    valid=all(k in t for k in ['userSeconds','systemSeconds'])
    dump(folder/'inner-resource.json',dict(**t,**status,
        cpuSeconds=t['userSeconds']+t['systemSeconds'] if valid else None,
        peakChildRssMiB=status.get('peakProcessRssMiB'),
        timingParseComplete=valid,rawTimingFile='inner-time.json',
        scope='bash time solver tree CPU; container process RSS summed every 0.5 s; malformed timing component retained as unavailable'))
    mesh_name=folder.name.split('-1e-')[0]
    mesh=BASE/'meshes'/f'{mesh_name}.solver.cm2d'
    if mesh.exists() and (folder/'solve.log').exists():
        annulus_cross.cases=lambda:{mesh_name:(mesh,Path(str(mesh).removesuffix('.solver.cm2d')+'.foam'))}
        result=collect_external(mesh_name,'simple',folder)
        result['solveSeconds']=t.get('wallSeconds');result['kernelTimeScope']='bash time around timeout/simpleFoam, excluding checkMesh';dump(folder/'result.json',result)
        if status.get('exitCode') not in (None,0) and result.get('converged'):
            result.update(converged=False,status='not-converged',reason='process exit was not successful')
            dump(folder/'result.json',result)

def same_mesh():
    records=[]
    for n in (128,256):
        left=BASE/'meshes'/f'graded{n}';right=BASE/'meshes'/f'graded-indexed{n}'
        files=['.raw.cm2d','.solver.cm2d','.boundaries']
        if all(Path(str(p)+s).exists() for p in (left,right) for s in files):
            records.append(dict(n=n,files={s:dict(baselineSha256=sha(Path(str(left)+s)),
                indexedSha256=sha(Path(str(right)+s)),identical=sha(Path(str(left)+s))==sha(Path(str(right)+s))) for s in files}))
    dump(BASE/'indexed-equivalence.json',records)
    pairs=[('shape160-coupled-adaptive-ic0-1e-8','shape160-coupled-adaptive-aggregation-1e-8-connected'),
           ('graded128-simple-strict-ic0-1e-6','graded128-coupled-adaptive-ic0-1e-8')]
    comparisons=[]
    for a,b in pairs:
        files=[BASE/'solutions'/f'{name}.cells.csv' for name in (a,b)]
        if not all(p.exists() for p in files):continue
        tables=[]
        for file in files:
            with file.open(encoding='utf-8-sig') as f:tables.append(list(csv.DictReader(f)))
        if len(tables[0])!=len(tables[1]):raise ValueError('same-mesh field count mismatch')
        comparisons.append(dict(first=a,second=b,cells=len(tables[0]),
            velocityComponentMaxOverUref=max(abs(float(x[k])-float(y[k]))/.5 for x,y in zip(*tables) for k in ['u','v']),
            pressureMaxDifference=max(abs(float(x['p'])-float(y['p'])) for x,y in zip(*tables)),
            scope='Measured same-cell iterative result differences, not new physical accuracy thresholds'))
    dump(BASE/'same-mesh-fields.json',comparisons)

if __name__=='__main__':
    for folder in (BASE/'fluent').glob('*'):
        if (folder/'result.json').exists():fluent_fields(folder)
    for folder in (BASE/'openfoam').glob('*'):
        if folder.is_dir():foam_resources(folder)
    same_mesh()
    print('postprocessing finished')
