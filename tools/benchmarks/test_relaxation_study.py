"""Research bookkeeping checks; actual PDE verification stays in native tests."""
import json
from pathlib import Path
import tempfile
import unittest

import relaxation_study as study


class StudyBookkeepingTests(unittest.TestCase):
    def test_failed_or_different_solution_cannot_win(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = dict(training=['a', 'b'])
            index = 0
            for case in ['a', 'b', 'heldout']:
                for pair, status, seconds in [([.6, .25], 'eligible', 10),
                                               ([1., 1.], 'different-solution', .1),
                                               ([.2, .1], 'time-limit', .01)]:
                    study.write(root / 'runs/grid' / str(index) / 'record.json',
                                dict(case=case, batch='grid', method='simple', pair=pair,
                                     status=status, wallSeconds=seconds, reference=False))
                    index += 1
            # A successful holdout cannot rescue a missing/failed training pair.
            study.write(root / 'runs/grid/extra/record.json',
                        dict(case='heldout', batch='grid', method='simple', pair=[.8, .8],
                             status='eligible', wallSeconds=.001, reference=False))
            for case in ['a', 'b']:
                study.write(root / f'runs/grid/warm-{case}/record.json',
                            dict(case=case, batch='grid', method='simple', pair=[.9,.9],
                                 status='eligible', wallSeconds=.001, reference=False,
                                 initialHashes={'cells':'warm-start-input'}))
            result = study.select(root, manifest, ['grid'], root/'selection.json')
            self.assertEqual([list(r['pair']) for r in result['rankings']['simple']], [[.6, .25]])
            # Reproducing a frozen selection must compare equal after JSON I/O.
            self.assertEqual(study.select(root,manifest,['grid'],root/'selection.json'),result)

    def test_pressure_gauge_and_face_identity(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); a = root/'a'; b = root/'b'
            for prefix, offset in [(a, 0), (b, 5)]:
                Path(str(prefix)+'.cells.csv').write_text(
                    'cell,x,y,area,u,v,p\n' + f'0,0,0,1,1,0,{offset}\n1,1,0,3,2,0,{offset+2}\n')
                Path(str(prefix)+'.faces.csv').write_text('face,owner,neighbour,flux\n0,0,1,1\n1,1,-1,0\n')
            case = dict(closed=True, speed=1, nu=.1, length=1)
            self.assertEqual(study.agreement(a,b,case)['pressure'], 0)
            self.assertGreater(study.agreement(a,b,{**case, 'closed': False})['pressure'], 4)
            Path(str(b)+'.faces.csv').write_text('face,owner,neighbour,flux\n1,1,-1,0\n0,0,1,1\n')
            with self.assertRaisesRegex(ValueError, 'Face identity'):
                study.agreement(a,b,case)

    def test_caps_and_strict_certification(self):
        valid = dict(converged=True, strictLinearFinal=True)
        errors = dict(finite=True, velocity=0, pressure=0, flux=0)
        self.assertEqual(study.classify(2,False,dict(converged=False),None,1e-4,''), 'iteration-limit')
        self.assertEqual(study.classify(-9,True,valid,errors,1e-4,''), 'time-limit')
        self.assertEqual(study.classify(0,False,{**valid, 'strictLinearFinal':False},errors,1e-4,''), 'uncertified')
        self.assertEqual(study.classify(0,False,valid,{**errors, 'flux':.1},1e-4,''), 'different-solution')

    def test_resume_cannot_change_experiment(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'spec.json'
            study.freeze(path, dict(pair=[.6,.25], seconds=60))
            study.freeze(path, dict(pair=[.6,.25], seconds=60))
            with self.assertRaisesRegex(ValueError, 'Immutable experiment changed'):
                study.freeze(path, dict(pair=[.6,.5], seconds=60))
            self.assertEqual(json.loads(path.read_text())['pair'], [.6,.25])


if __name__ == '__main__':
    unittest.main()
