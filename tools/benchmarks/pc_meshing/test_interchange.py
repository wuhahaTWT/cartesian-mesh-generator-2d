"""Minimal geometry/interchange regressions; no independent FV equations."""
from pathlib import Path
import csv,json,math,re,sys,tempfile,unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parent))
import h5py,numpy as np
from fluent_topology import case_cell_map
from fluent_reference import export_mesh
from postprocess import fluent_fields
import native_mesh as nm

def fixture(path):
    vertices=((1.,1.),(2.,1.),(4.,1.),(1.,2.),(2.,2.),(4.,2.))
    # Deliberately reversed stored edge 2: the consumer must recover owner orientation.
    edges=(nm.Edge(0,0,1,0,-1,2),nm.Edge(1,1,4,0,1,0),nm.Edge(2,3,4,0,-1,2),
           nm.Edge(3,3,0,0,-1,2),nm.Edge(4,1,2,1,-1,2),nm.Edge(5,2,5,1,-1,2),
           nm.Edge(6,5,4,1,-1,2))
    cells=(nm.Cell(0,1.,(0,1,4,3),(0,1,2,3)),nm.Cell(1,2.,(1,2,5,4),(4,5,6,1)))
    return nm.Mesh(path,vertices,edges,cells,())

def write_case(path,mesh):
    ordered=[e for e in mesh.edges if e.neighbour>=0]+[e for e in mesh.edges if e.neighbour<0]
    with h5py.File(path,'w') as f:
        b=f.create_group('meshes/1')
        b.create_dataset('nodes/coords/1',data=np.asarray(mesh.vertices))
        b.create_dataset('faces/nodes/1/nnodes',data=np.full(len(ordered),2))
        b.create_dataset('faces/nodes/1/nodes',data=np.asarray([[e.v0+1,e.v1+1] for e in ordered]).ravel())
        c0=b.create_dataset('faces/c0/1',data=np.asarray([2 if e.owner==0 else 1 for e in ordered]))
        c0.attrs['minId']=[1]
        c1=b.create_dataset('faces/c1/1',data=np.asarray([1 for e in ordered if e.neighbour>=0]))
        c1.attrs['minId']=[1]

class InterchangeTests(unittest.TestCase):
    def test_fluent_2d_positive_cell_orientation(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);source=root/'unit.cm2d';source.write_text('synthetic two-cell fixture')
            mesh=fixture(source)
            with patch('fluent_reference.nm.read_cm2d',return_value=mesh):
                export_mesh(source,root/'unit.msh')
            text=(root/'unit.msh').read_text(encoding='utf-8')
            blocks=re.findall(r'\(13 \([^\n]+\)\(\n(.*?)\)\)',text,re.S)
            signed=[0.,0.]
            for block in blocks:
                for line in block.strip().splitlines():
                    a,b,c0,c1=[int(x,16) for x in line.split()]
                    x0,y0=mesh.vertices[a-1];x1,y1=mesh.vertices[b-1]
                    contribution=(x0*y1-y0*x1)/2
                    if c0:signed[c0-1]+=contribution
                    if c1:signed[c1-1]-=contribution
            self.assertEqual(signed,[1.,2.])

    def test_parallel_fluent_export_uses_exact_case_topology_not_export_ids(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);source=root/'unit.cm2d';source.write_text('synthetic two-cell fixture')
            mesh=fixture(source);g=nm.measure(mesh)
            (root/'result.json').write_text(json.dumps(dict(status='finished',mesh=str(source))))
            with (root/'cells.csv').open('w',newline='') as f:
                w=csv.writer(f);w.writerow(['cellnumber','x-coordinate','y-coordinate','x-velocity','y-velocity','pressure'])
                w.writerow([1,*g.centroids[0],.5,.6,5.7])
                w.writerow([2,*g.centroids[1],.2,.3,5.4])
            write_case(root/'accepted.cas.h5',mesh)
            with h5py.File(root/'accepted.dat.h5','w') as f:
                f.create_dataset('settings/Case File',data=np.array([b'accepted.cas.h5']))
                for name,values in [('SV_U',[.2,.5]),('SV_V',[.3,.6]),('SV_P',[.4,.7])]:
                    d=f.create_dataset('results/1/phase-1/cells/'+name+'/1',data=values)
                    d.attrs['minId']=[1];d.attrs['maxId']=[2]
            with patch('postprocess.nm.read_cm2d',return_value=mesh):fluent_fields(root)
            record=json.loads((root/'result.json').read_text(encoding='utf-8'))
            self.assertTrue(record['cellMapping']['exactNodeAndCellIncidenceBijection'])
            self.assertFalse(record['cellMapping']['caseIdsPreserved'])
            with (root/'comparison.cells.csv').open() as f:rows=list(csv.DictReader(f))
            self.assertEqual([float(r['area']) for r in rows],[1.,2.])
            self.assertEqual([float(r['u']) for r in rows],[.5,.2])

    def test_changed_wall_coordinate_cannot_match(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);mesh=fixture(root/'unit.cm2d');write_case(root/'unit.cas.h5',mesh)
            with h5py.File(root/'unit.cas.h5','a') as f:
                f['meshes/1/nodes/coords/1'][0,0]=1.01
            with self.assertRaisesRegex(ValueError,'changed a node'):
                case_cell_map(mesh,root/'unit.cas.h5')

if __name__=='__main__':unittest.main()
