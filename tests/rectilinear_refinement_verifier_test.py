import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools/verification'))
import verify_rectilinear_refinement as v


def grid(x,y):
    return [[(a,c),(b,c),(b,d),(a,d)] for c,d in zip(y,y[1:]) for a,b in zip(x,x[1:])]


class NestedRefinement(unittest.TestCase):
    def test_partition_and_preserved_band(self):
        a=grid([-.5,0,1],[0,.125,.25,1])
        b=grid([-.5,0,1],[0,.125,.1875,.25,.625,1])
        # Cell order and polygon starting vertex are irrelevant to geometry.
        b=[p[2:]+p[:2] for p in reversed(b)]
        r=v.check_refinement(a,b,1)
        self.assertEqual((r['originalCells'],r['refinedCells'],r['preservedCells']),(6,10,2))
        self.assertEqual(r['domainArea'],1.5)
        self.assertTrue(r['areaIdenticalAsExactRational'])

    def test_invalid_refinement_rejected(self):
        a=grid([0,.5,1],[0,.125,.25,1])
        good=grid([0,.5,1],[0,.125,.1875,.25,.625,1])
        for bad in [good[:-1],good+[good[0]],
                    grid([0,.6,1],[0,.125,.1875,.25,.625,1]),
                    grid([0,.5,1],[0,.0625,.125,.1875,.25,.625,1]),
                    grid([0,.5,1],[0,.125,.19,.25,.625,1]),
                    grid([0,.5,1],[0,.125,.1875,.25,.625,2])]:
            with self.assertRaises(ValueError):v.check_refinement(a,bad,1)
        for count in (0,3,True,1.5):
            with self.assertRaises(ValueError):v.check_refinement(a,good,count)
        with self.assertRaises(ValueError):v.tensor_axes([[(0,0),(1,0),(1,1),(.1,1)]])
        with self.assertRaises(ValueError):v.tensor_axes([[(0,0),(1,1),(1,0),(0,1)]])
        with self.assertRaises(ValueError):v.tensor_axes([[(0,0),(1,0),(1,1),(0,float('nan'))]])


if __name__=='__main__':unittest.main()
