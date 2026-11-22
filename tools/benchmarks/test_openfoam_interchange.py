"""Minimal regressions for real foam-extend interchange failures; no PDE audit."""
import tempfile
import unittest
from pathlib import Path
import openfoam as reference

foam=reference.foam

class Interchange(unittest.TestCase):
    def test_signed_legacy_neighbours_with_footer(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'neighbour'
            p.write_text(foam.header('labelList','constant/polyMesh','neighbour')+
                         '3\n(\n7\n-1\n-1\n)\n// legacy footer\n')
            self.assertEqual(foam.read_labels(p),[7,-1,-1])

    def fixture(self,root,full):
        source=root/'source/constant/polyMesh';source.mkdir(parents=True)
        # A small extruded annular sector exercises physical and empty patches.
        # This is a file-interchange fixture, never a flow-accuracy benchmark.
        xy=[(.49361353808843056,.07942309080391574),(.4969160653177267,.0551069538821942),
            (.9938321306354534,.1102139077643884),(.9872270761768611,.15884618160783148)]
        points=[(x,y,z) for z in (-.05,.05) for x,y in xy]
        faces=[(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7),(3,2,1,0),(4,5,6,7)]
        def write(name,cls,records):
            (source/name).write_text(foam.header(cls,'constant/polyMesh',name)+str(len(records))+'\n(\n'+'\n'.join(records)+'\n)\n// legacy footer\n')
        write('points','vectorField',['('+' '.join(f'{x:.17g}' for x in p)+')' for p in points])
        write('faces','faceList',[f'{len(f)}('+' '.join(map(str,f))+')' for f in faces])
        write('owner','labelList',['0']*6)
        write('neighbour','labelList',['-1']*6 if full else [])
        write('boundary','polyBoundaryMesh',[
            'rotor { type patch; nFaces 1; startFace 0; }',
            'housing { type patch; nFaces 3; startFace 1; }',
            'caps { type empty; nFaces 2; startFace 4; }'])
        return source.parents[1]

    def check_prepared(self,full):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=self.fixture(root,full);mesh=root/'input.cm2d';mesh.write_text('hash fixture')
            output=root/'case'
            reference.prepare('annulus',source,mesh,output,'upwind',20)
            patches=foam.read_boundary(output/'constant/polyMesh/boundary')
            self.assertEqual(sum(p['nFaces'] for p in patches if p['type']=='empty'),2)
            self.assertGreater(sum(p['nFaces'] for p in patches if p['name']=='rotor'),0)
            self.assertGreater(sum(p['nFaces'] for p in patches if p['name']=='housing'),0)
            self.assertEqual(sum(p['nFaces'] for p in patches if p['type']!='empty'),4)
            self.assertEqual(foam.read_labels(output/'constant/polyMesh/neighbour'),[])
            self.assertIn('rotor\n{ type fixedValue;', (output/'0/U').read_text())

    def test_full_neighbour_normalization(self):self.check_prepared(True)
    def test_abbreviated_neighbour_preserves_patches(self):self.check_prepared(False)

if __name__=='__main__':unittest.main()
