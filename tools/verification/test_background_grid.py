#!/usr/bin/env python3
import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from check_background_grid import audit


def main():
    cli = str(Path(sys.argv[1]).resolve())
    with tempfile.TemporaryDirectory(prefix='cartmesh-background-') as folder:
        root = Path(folder)
        boundary = root/'body.xy'
        # Nested loops and a sloping boundary exercise parity and crossing cells.
        boundary.write_text('0 0\n2 0\n2 2\n0 2\n\n0.6 0.6\n1.4 0.6\n1.0 1.4\n')
        for mode in ('uniform', 'adaptive'):
            prefix = root/mode
            args = [cli,str(boundary),str(prefix),'5','0.5','0.1','exterior','-','2','--background-grid',mode]
            subprocess.run(args, check=True, capture_output=True, timeout=20)
            path = Path(str(prefix)+'.background.json')
            result = audit(path)
            assert all(n > 0 for n in result['counts'])
            original = path.read_bytes()
            vtk = Path(str(prefix)+'.background.vtk').read_bytes()
            # Check the second export independently, including cell-associated
            # classes and levels rather than only trusting JSON or exit status.
            tokens = vtk.decode().split()
            cursor = tokens.index('POINTS')
            n = len(json.loads(original)['cells'])
            assert int(tokens[cursor+1]) == 4*n and tokens[cursor+2] == 'double'
            points = list(map(float, tokens[cursor+3:cursor+3+12*n]))
            cursor += 3+12*n
            assert tokens[cursor:cursor+3] == ['CELLS',str(n),str(5*n)]
            cursor += 3
            for i, cell in enumerate(json.loads(original)['cells']):
                assert list(map(int,tokens[cursor:cursor+5])) == [4,4*i,4*i+1,4*i+2,4*i+3]
                cursor += 5
                x0,y0,x1,y1 = cell['bounds']
                assert points[12*i:12*i+12] == [x0,y0,0,x1,y0,0,x1,y1,0,x0,y1,0]
            assert tokens[cursor:cursor+2] == ['CELL_TYPES',str(n)]
            cursor += 2
            assert tokens[cursor:cursor+n] == ['9']*n
            cursor += n
            assert tokens[cursor:cursor+2] == ['CELL_DATA',str(n)]
            cursor += 2
            for field in ('classification', 'level'):
                assert tokens[cursor:cursor+6] == ['SCALARS',field,'int','1','LOOKUP_TABLE','default']
                cursor += 6
                assert list(map(int,tokens[cursor:cursor+n])) == [c[field] for c in json.loads(original)['cells']]
                cursor += n
            assert cursor == len(tokens)
            subprocess.run(args, check=True, capture_output=True, timeout=20)
            assert path.read_bytes() == original
            assert Path(str(prefix)+'.background.vtk').read_bytes() == vtk
            assert not Path(str(prefix)+'.cm2d').exists()
            data = json.loads(original)
            for mutation in ('classification', 'missing', 'duplicate', 'overlap'):
                bad = copy.deepcopy(data)
                if mutation == 'classification': bad['cells'][0]['classification'] = (bad['cells'][0]['classification']+1)%3
                elif mutation == 'missing': bad['cells'].pop()
                elif mutation == 'duplicate': bad['cells'].append(copy.deepcopy(bad['cells'][0]))
                else:
                    bad['cells'][0].update(level=0, ix=0, iy=0, bounds=bad['domain'])
                bad_path = root/'bad.json'
                bad_path.write_text(json.dumps(bad))
                try: audit(bad_path)
                except AssertionError: pass
                else: raise AssertionError(f'accepted corrupt {mutation}')
        for extra in (['--background-grid','invalid'], ['--background-grid','uniform','--background-grid','uniform']):
            run = subprocess.run([cli,str(boundary),str(root/'bad'),'5','0.5','0.1','exterior','-','2',*extra], capture_output=True, timeout=20)
            assert run.returncode != 0
        for region, output, depth in (('interior','-',5), ('exterior',str(root/'foam'),5), ('exterior','-',11)):
            run = subprocess.run([cli,str(boundary),str(root/'bad'),str(depth),'0.5','0.1',region,output,'2','--background-grid','uniform'], capture_output=True, timeout=20)
            assert run.returncode != 0
        boundary.write_text('0 0\n1 1\n0 1\n1 0\n')
        run = subprocess.run([cli,str(boundary),str(root/'bad'),'5','--background-grid','adaptive'], capture_output=True, timeout=20)
        assert run.returncode != 0
    print('Background uniform/adaptive coverage, classification, balance, deterministic export and refusal checks passed')


if __name__ == '__main__': main()
