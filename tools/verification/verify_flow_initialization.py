"""Independent initial-field reconstruction from CM2D and a retained checkpoint."""
import math
import shlex
from pathlib import Path
import verify_native_flow as native


def fail(message):
    raise native.VerificationError('initial vortex: '+message)


def audit(prefix, mesh, measured, summary, history, cells):
    spec = summary.get('initialVortex')
    if spec is None:
        if Path(str(prefix)+'.initial.checkpoint').exists():
            fail('initial checkpoint has no defining metadata')
        return None
    if not isinstance(spec, dict) or set(spec) != {'definition','centre','radius','peakSpeed','checkpointSuffix'}:
        fail('malformed definition')
    if spec['definition'] != 'compact-cubic-v1' or spec['checkpointSuffix'] != '.initial.checkpoint':
        fail('unsupported definition/checkpoint suffix')
    if summary.get('case') not in ('external','channel','duct','cavity','custom') or abs(history[0]['time']-history[0]['dt'])>1e-14:
        fail('requires a fresh physical transient calculation')
    if not isinstance(spec['centre'], list) or len(spec['centre']) != 2:
        fail('invalid centre')
    cx, cy = (native.finite(v,'initial centre') for v in spec['centre'])
    radius = native.finite(spec['radius'],'initial radius')
    speed = native.finite(spec['peakSpeed'],'initial peak speed')
    if radius <= 0: fail('radius must be positive')
    geometry = native.face_geometry(mesh, measured)
    ends = []
    winding, clearance = 0, math.inf
    for edge in mesh.edges:
        cell = mesh.cells[edge.owner]
        j = cell.edges.index(edge.id)
        a, b = mesh.vertices[cell.vertices[j]], mesh.vertices[cell.vertices[(j+1)%len(cell.vertices)]]
        ends.append((a,b))
        if edge.neighbour >= 0: continue
        dx,dy = b[0]-a[0],b[1]-a[1]
        qx,qy = cx-a[0],cy-a[1]
        t = max(0,min(1,(qx*dx+qy*dy)/(dx*dx+dy*dy)))
        clearance = min(clearance,math.hypot(qx-t*dx,qy-t*dy))
        cross = dx*qy-dy*qx
        winding += int(a[1]<=cy<b[1] and cross>0)-int(b[1]<=cy<a[1] and cross<0)
    if winding != 1 or clearance <= radius: fail('support intersects a boundary or is outside fluid')
    path = Path(str(prefix)+'.initial.checkpoint')
    rows = iter(shlex.split(line) for line in path.read_text().splitlines() if line.strip())
    def row(label, length=None):
        value = next(rows, [])
        if not value or value[0] != label or (length is not None and len(value)!=length):
            fail('malformed '+label)
        return value[1:]
    def numbers(values): return [native.finite(v,'initial checkpoint value') for v in values]
    def equal(values, expected, label):
        actual = numbers(values)
        if len(actual)!=len(expected) or any(not native.close(a,b,1e-12,1e-10) for a,b in zip(actual,expected)):
            fail('checkpoint '+label+' differs from reconstructed input')
    version = row('CARTMESH2D_FLOW_CHECKPOINT',2)[0]
    if version not in ('2','3','4') or (version=='4') != (summary['case']=='custom'): fail('unsupported checkpoint version')
    if row('DISCRETIZATION',2) != ['Euler-RC-v2']: fail('unsupported discretization')
    config = row('CONFIG',8)
    if (config[0],config[3],config[4],config[6]) != (summary['case'],summary['convection'],summary['viscousStress'],summary['outletBackflow']):
        fail('checkpoint physics configuration differs')
    equal([config[1],config[2],config[5]],[summary['nu'],summary['speed'],0],'configuration')
    if version in ('3','4'):
        viscosity = row('FACE_VISCOSITY')
        n = native.integer(viscosity[0],'viscosity count')
        if len(viscosity)!=n+1 or n not in (0,len(mesh.edges)): fail('bad viscosity count')
        if n:
            source = native.load_csv(Path(summary['viscosityFile']),('face','viscosity'))
            values = {native.integer(r['face'],'viscosity face'):native.finite(r['viscosity'],'viscosity') for r in source}
            if set(values)!=set(range(n)) or min(values.values())<=0: fail('invalid viscosity input')
            equal(viscosity[1:],[values[i] for i in range(n)],'face viscosity')
    initial_pressure = 0
    if version=='4':
        boundaries = native.audit_explicit_boundaries(prefix,mesh,measured,summary)
        if row('BOUNDARIES',2)!=[str(len(boundaries))]: fail('boundary count differs')
        for boundary in boundaries:
            fields = row('BOUNDARY',7)
            if fields[:3]!=[str(boundary['face']),boundary['type'],boundary['name']]: fail('boundary identity differs')
            equal(fields[3:],[boundary[k] for k in ('u','v','p')],'boundary values')
        outlets = [(math.hypot(*geometry[b['face']].area_vector),b['p']) for b in boundaries if b['type']=='pressure-outlet']
        if outlets: initial_pressure = math.fsum(w*p for w,p in outlets)/math.fsum(w for w,p in outlets)
    if row('CELLS',2)!=[str(len(mesh.cells))]: fail('cell count differs')
    for cell in mesh.cells:
        fields = row('CELL',6+len(cell.edges))
        if fields[0]!=str(cell.id) or fields[4:]!=[str(len(cell.edges)),*map(str,cell.edges)]: fail('cell topology differs')
        equal(fields[1:4],[*measured.centroids[cell.id],measured.areas[cell.id]],'cell geometry')
    if row('FACES',2)!=[str(len(mesh.edges))]: fail('face count differs')
    for edge,g in zip(mesh.edges,geometry):
        fields = row('FACE',13)
        if fields[:4]!=[str(edge.id),str(edge.owner),str(edge.neighbour) if edge.neighbour>=0 else '-',str(edge.patch)]: fail('face topology differs')
        equal(fields[4:],[*g.centre,*g.area_vector,*g.correction,g.transmissibility,g.neighbour_weight],'face geometry')
    if row('TIME',2)!=['0']: fail('initial physical time is not zero')
    state = {}
    for key,n in (('U',len(mesh.cells)),('V',len(mesh.cells)),('P',len(mesh.cells)),('FLUX',len(mesh.edges))):
        fields = row(key,n+2)
        if fields[0]!=str(n): fail('state count differs')
        state[key] = numbers(fields[1:])
    row('END',1)
    if next(rows,None) is not None: fail('trailing checkpoint data')
    def sample(point):
        x,y = (point[0]-cx)/radius,(point[1]-cy)/radius
        q = x*x+y*y
        if q>=1: return 0.,0.,0.
        s = 1-q
        a = speed*25*math.sqrt(5)/16
        return -a*y*s*s,a*x*s*s,a*radius*s*s*s/6
    expected = [sample(point) for point in measured.centroids]
    equal(state['U'],[v[0] for v in expected],'initial U')
    equal(state['V'],[v[1] for v in expected],'initial V')
    equal(state['P'],[initial_pressure]*len(mesh.cells),'initial P')
    flux = [sample(b)[2]-sample(a)[2] if e.neighbour>=0 else 0. for e,(a,b) in zip(mesh.edges,ends)]
    equal(state['FLUX'],flux,'initial FLUX')
    if speed!=0 and (not any(state['U']+state['V']) or not any(flux)):
        fail('nonzero seed is unresolved by cells or faces')
    if len(history)==1:
        equal([c['previousU'] for c in cells],state['U'],'first-step previous U')
        equal([c['previousV'] for c in cells],state['V'],'first-step previous V')
    balances = [math.fsum(state['FLUX'][f]*(1 if mesh.edges[f].owner==c.id else -1) for f in c.edges) for c in mesh.cells]
    scale = max(abs(speed)*radius,1e-300)
    if max(map(abs,balances))>1e-12*scale: fail('initial face flux is not locally conservative')
    return dict(valid=True,definition=spec,checkpointSha256=native.sha256_file(path),
        minimumBoundaryClearance=clearance,maximumCellFluxImbalance=max(map(abs,balances)),
        sampledPeakSpeed=max(math.hypot(u,v) for u,v in zip(state['U'],state['V'])))
