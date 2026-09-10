import math
from pathlib import Path
import sys
import unittest
import tempfile
import subprocess

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools' / 'verification'))
from check_face_planes import classify, face_geometry, run_case

def prism(poly):
    n=len(poly); pts=[(x,y,-.01) for x,y in poly]+[(x,y,.01) for x,y in poly]
    fs=[list(range(n-1,-1,-1)),list(range(n,2*n))]
    for i in range(n): fs.append([i,(i+1)%n,(i+1)%n+n,i+n])
    return pts,fs,[0]*len(fs)

def edge_balance(faces):
    edges={}
    for f in faces:
        for a,b in zip(f,f[1:]+f[:1]): edges.setdefault(tuple(sorted((a,b))),[]).append((a,b))
    return all(len(v)==2 and v[0]==(v[1][1],v[1][0]) for v in edges.values())

class FacePlaneAuditTest(unittest.TestCase):
    def test_convex_coplanar_and_concave_classification(self):
        p,f,o=prism([(-1,-1),(1,-1),(1,1),(-1,1)])
        self.assertTrue(edge_balance(f)); self.assertEqual(classify([list(range(len(f)))],p,f,o)['safe'],{0})
        p,f,o=prism([(-1,-1),(1,-1),(1,0),(1,1),(-1,1)])
        self.assertTrue(edge_balance(f)); self.assertEqual(classify([list(range(len(f)))],p,f,o)['near_coplanar_only'],{0})
        p,f,o=prism([(-1,-1),(1,-1),(1,0),(0,0),(0,1),(-1,1)])
        self.assertTrue(edge_balance(f)); self.assertEqual(classify([list(range(len(f)))],p,f,o)['positive_side'],{0})

    def test_two_cell_owner_neighbour_flip(self):
        p=[(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1),(2,0,0),(2,1,0),(2,0,1),(2,1,1)]
        # Internal face is first, outward from owner 0 toward neighbour 1.
        f=[[1,2,6,5], [0,3,2,1], [4,5,6,7], [0,4,7,3],
           [0,1,5,4], [3,7,6,2], [1,2,9,8], [5,10,11,6],
           [8,9,11,10], [1,8,10,5], [2,6,11,9]]
        o=[0]*6+[1]*5
        cells=[list(range(6)),[0,6,7,8,9,10]]
        self.assertEqual((len(p),len(f)),(12,11))
        for ci,fs in enumerate(cells):
            self.assertEqual(len(fs),6)
            self.assertTrue(edge_balance([f[i] if o[i]==ci else f[i][::-1] for i in fs]))
        self.assertEqual(classify(cells,p,f,o)['safe'],{0,1})
        f[0]=list(reversed(f[0])); o[0]=1
        self.assertEqual(classify(cells,p,f,o)['safe'],{0,1})
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            (root/'points').write_text('12\n(\n'+'\n'.join('(%g %g %g)'%x for x in p)+'\n)\n')
            (root/'faces').write_text('11\n(\n'+'\n'.join('%d(%s)'%(len(x),' '.join(map(str,x))) for x in f)+'\n)\n')
            (root/'owner').write_text('11\n('+' '.join(map(str,o))+')\n')
            (root/'neighbour').write_text('1\n(0)\n')
            self.assertEqual(run_case(root)['safe'],2)

    def test_reject_bad_indices_and_nonfinite(self):
        p,f,o=prism([(-1,-1),(1,-1),(1,1),(-1,1)])
        with self.assertRaises(ValueError): face_geometry([999],p)
        with self.assertRaises(ValueError): face_geometry([0,1,2],[(0,0,0),(math.inf,0,0),(0,1,0)])
        with self.assertRaises(ValueError): classify([[0]],p,f[:1], [0])
        with self.assertRaises(ValueError): face_geometry([],p)
        with self.assertRaises(ValueError): face_geometry([0,1,1],p)
        with self.assertRaises(ValueError): classify([[0,1,2,2]],p,f,o)
        with self.assertRaises(ValueError): classify([[0,1,2,3]],p,f,o)
        with self.assertRaises(ValueError): classify([list(range(6))],p,f,[1]*6)

    def test_serialized_input_validation_and_expected_set_bounds(self):
        p,f,o=prism([(-1,-1),(1,-1),(1,1),(-1,1)])
        def vecs(): return '8\n(\n'+'\n'.join('(%g %g %g)'%x for x in p)+'\n)\n'
        def flist(): return '6\n(\n'+'\n'.join('%d(%s)'%(len(x),' '.join(map(str,x))) for x in f)+'\n)\n'
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); (root/'points').write_text(vecs()); (root/'faces').write_text(flist()); (root/'owner').write_text('6\n(0 0 0 0 0 0)\n'); (root/'neighbour').write_text('0\n(\n)\n')
            self.assertEqual(run_case(root)['safe'],1)
            wrong=root/'wrong-set'; wrong.write_text('1\n(0)\n')
            self.assertFalse(run_case(root,wrong)['expected_set_matches'])
            proc=subprocess.run([sys.executable,str(Path(__file__).resolve().parents[1]/
                'tools/verification/check_face_planes.py'),str(root),'--expected-set',str(wrong)],
                capture_output=True,text=True)
            self.assertEqual(proc.returncode,1)
            exp=root/'set'; exp.write_text('1\n(2)\n')
            with self.assertRaises(ValueError): run_case(root,exp)
            (root/'points').write_text(vecs().replace('8\n','7\n',1))
            with self.assertRaises(ValueError): run_case(root)
            (root/'points').write_text(vecs().replace('-1 -1 -0.01','nan -1 -0.01',1))
            with self.assertRaises(ValueError): run_case(root)
            (root/'points').write_text(vecs())
            (root/'owner').write_text('6\n(1 0 0 0 0 0)\n')
            with self.assertRaises(ValueError): run_case(root)
            (root/'owner').write_text('6\n(0 0 0 0 0 0)\n'); (root/'neighbour').write_text('1\n(0)\n')
            with self.assertRaises(ValueError): run_case(root)
            (root/'neighbour').write_text('1\n(-1)\n')
            with self.assertRaises(ValueError): run_case(root)

    def test_serialized_garbage_and_inconsistent_connectivity(self):
        p,f,o=prism([(-1,-1),(1,-1),(1,1),(-1,1)])
        baseline={
            'points':'8\n(\n'+'\n'.join('(%g %g %g)'%x for x in p)+'\n)\n',
            'faces':'6\n(\n'+'\n'.join('%d(%s)'%(len(x),' '.join(map(str,x))) for x in f)+'\n)\n',
            'owner':'6\n(0 0 0 0 0 0)\n', 'neighbour':'0\n()\n'}
        corrupt=[('owner','5\n(0 0 0 0 0)\n'),
                 ('owner','6\n(0oops 0 0 0 0 0)\n'),
                 ('owner','6\n(0 0 0 0 0 999999999)\n'),
                 ('neighbour','1\n(0)\n'), ('neighbour','1\n(-1)\n'),
                 ('neighbour','2\n(1)\n'),
                 ('faces',baseline['faces'].replace('4(', '5(',1)),
                 ('faces',baseline['faces'].replace('4(', '4(junk ',1)),
                 ('faces',baseline['faces'].replace('3 2 1 0','3 2 1 1',1)),
                 ('points',baseline['points'].replace('8\n(', '8\n(junk',1)),
                 ('points',baseline['points']+'junk'),
                 ('points',baseline['points'].replace('-1 -1 -0.01','nan -1 -0.01',1))]
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            for name,text in corrupt:
                with self.subTest(file=name,text=text):
                    for key,value in baseline.items(): (root/key).write_text(value)
                    (root/name).write_text(text)
                    with self.assertRaises(ValueError): run_case(root)

if __name__=='__main__': unittest.main()
