#!/usr/bin/env python3
"""Independent analytic interpolation/scaling and rejection tests."""
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools/verification'))
import analyze_sst_flatplate as analysis


class FlatPlateAnalysis(unittest.TestCase):
    def test_affine_profile_and_wall_units(self):
        x,y=[.25,.75],[.01,.04,.15]
        grid={(i,j):{'u':2*xx+3*yy,'v':0.,'p':1.,'k':.01,'omega':2.,'nuT':.001}
              for i,xx in enumerate(x) for j,yy in enumerate(y)}
        wall=[{'x':.25,'kinematicShear':.04},{'x':.75,'kinematicShear':.08}]
        result=analysis.station_profile(x,y,grid,wall,.5,.01,2.)
        self.assertAlmostEqual(result['Cf'],.03)
        self.assertAlmostEqual(result['uTau'],.06**.5)
        for yy,row in zip(y,result['samples']):
            self.assertAlmostEqual(row['u'],1+3*yy)
            self.assertAlmostEqual(row['uOverU'],(1+3*yy)/2)
            self.assertAlmostEqual(row['uPlus'],(1+3*yy)/.06**.5)
            self.assertAlmostEqual(row['yPlus'],yy*.06**.5/.01)
        self.assertIsNone(result['delta99FirstCrossing'])
        with self.assertRaises(ValueError):analysis.station_profile(x,y,grid,wall,.1,.01,2.)
        for sample in wall:sample['kinematicShear']=-.1
        with self.assertRaises(ValueError):analysis.station_profile(x,y,grid,wall,.5,.01,2.)

    def test_no_extrapolation_or_invalid_coordinates(self):
        self.assertEqual(analysis.interpolate([0,2],[1,5],0),1)
        self.assertEqual(analysis.interpolate([0,2],[1,5],2),5)
        self.assertEqual(analysis.interpolate([0,2],[1,5],.5),2)
        for x,y,t in [([0,0],[1,2],0),([0,2],[1,float('nan')],1),([0,2],[1,2],3),([0],[1],0)]:
            with self.assertRaises(ValueError):analysis.interpolate(x,y,t)
        # Overshoot/nonmonotonic profiles use the FIRST 99% crossing; absence is explicit.
        self.assertAlmostEqual(analysis.first_crossing([0,1,2,3],[0,1.1,.8,1]),.9)
        self.assertIsNone(analysis.first_crossing([0,1,2],[0,.4,.8]))
        with self.assertRaises(ValueError):analysis.first_crossing([],[])
        with self.assertRaises(ValueError):analysis.first_crossing([0,1],[0,1],0)

    def test_actual_polygon_grid_required(self):
        vertices=[(0,0),(1,0),(0,1),(1,1),(0,2),(1,2)]
        cells=[SimpleNamespace(id=0,vertices=(0,1,3,2)),SimpleNamespace(id=1,vertices=(2,3,5,4))]
        mesh=SimpleNamespace(vertices=vertices,cells=cells)
        measured=SimpleNamespace(centroids=[(.5,.5),(.5,1.5)])
        rows=[dict(cell=i,u=i,v=0,p=0,k=1,omega=1,nuT=1) for i in range(2)]
        xs,ys,grid=analysis.tensor_field(mesh,measured,rows)
        self.assertEqual((xs,ys),([.5],[.5,1.5]))
        self.assertEqual(grid[0,1]['u'],1)
        cells[1]=SimpleNamespace(id=1,vertices=(0,1,3,2))
        with self.assertRaises(ValueError):analysis.tensor_field(mesh,measured,rows)
        cells[1]=SimpleNamespace(id=1,vertices=(2,3,5))
        with self.assertRaises(ValueError):analysis.tensor_field(mesh,measured,rows)

    def test_reference_table_rejection(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'reference.dat'
            p.write_text('variables="x","u"\nzone,t="A"\n0 1\n2 5\n')
            self.assertEqual(analysis.read_zones(p,2),{'A':[[0,1],[2,5]]})
            for text in ['0 1\n2 5\n','zone,t="A"\n0 nan\n2 5\n',
                         'zone,t="A"\n0 1\n2 5\nzone,t="A"\n0 1\n2 5\n']:
                p.write_text(text)
                with self.assertRaises(ValueError):analysis.read_zones(p,2)

    def test_comparison_uses_common_physical_heights(self):
        common={'physicalInputs':{'nu':.01},'model':'test','leadingEdge':0,'topBoundary':'test','convection':'test'}
        def case(name,ys,offset,cf):
            return {**common,'name':name,'CdSkinFriction':cf,
                    'profiles':[{'x':1,'Cf':cf,'samples':[{'y':y,'uOverU':2*y+offset} for y in ys]}]}
        a=case('a',[.001,.02,.1],0,.01);b=case('b',[.002,.03,.2],.005,.0102)
        result=analysis.compare_cases(a,b)['stations'][0]
        self.assertAlmostEqual(result['signedRelativeCfChange'],.02)
        self.assertAlmostEqual(result['maxAbsoluteVelocityChangeOverU'],.005)
        self.assertEqual(result['profileCommonY'][0],.002)
        self.assertEqual(result['profileCommonY'][-1],.05)
        b['physicalInputs']={'nu':.02}
        with self.assertRaises(ValueError):analysis.compare_cases(a,b)


if __name__=='__main__':unittest.main()
