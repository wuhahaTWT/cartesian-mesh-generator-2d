#!/usr/bin/env python3
"""Run the packaged desktop, then independently read its exported ZIP and mesh."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile

from check_background_grid import audit as audit_background

ROOT = Path(__file__).resolve().parents[2]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--app', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--electron-arg', action='append', default=[])
    args = parser.parse_args()
    output = args.out.resolve()
    output.mkdir(parents=True, exist_ok=True)
    results = []
    # Exercise UTF-8 plus spaces through Node, native argv/filesystem and ZIP.
    workspace = output / '中文 路径'
    workspace.mkdir(exist_ok=True)
    temporary = workspace / 'temp'; temporary.mkdir(exist_ok=True)
    env = dict(os.environ, TEMP=str(temporary), TMP=str(temporary), TMPDIR=str(temporary))
    for name, image, method in [('png', 'raster-input-L.png', 'cutcell'),
                                ('jpg', 'raster-input-L.jpg', 'cutcell'),
                                ('hybrid', None, 'hybrid'),
                                ('background-uniform', None, 'background'),
                                ('background-adaptive', None, 'background')]:
        target = workspace / name; target.mkdir(exist_ok=True)
        screenshot = target / 'app.png'
        archive = target / 'result.zip'
        command = [str(args.app.resolve()), *args.electron_arg, '--smoke=circle',
                   '--method=' + method, '--out=' + str(target / 'cases'),
                   '--export=' + str(archive), '--shot=' + str(screenshot)]
        if method == 'background':
            command += ['--background-mode=' + name.split('-')[1],
                        '--background-level=7', '--background-minimum-level=4', '--background-padding=.5']
        if image:
            source = target / ('输入 图片' + Path(image).suffix)
            shutil.copyfile(ROOT / 'artifacts/current' / image, source)
            command += ['--image=' + str(source), '--image-width=200']
        with (target / 'app.log').open('w', encoding='utf-8') as log:
            completed = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=240)
        if completed.returncode: raise RuntimeError(f'{name}: App failed; see {target / "app.log"}')
        ui = json.loads(Path(str(screenshot) + '.json').read_text(encoding='utf-8'))
        assert ui['bundledChineseFontLoaded'], 'Chinese interface font did not load'
        unpacked = target / 'unpacked'
        with zipfile.ZipFile(archive) as z:
            assert z.testzip() is None
            z.extractall(unpacked)
        preview = next(unpacked.rglob('mesh-preview.png'))
        assert preview.read_bytes()[:8] == b'\x89PNG\r\n\x1a\n'
        assert next(unpacked.rglob('README_CN.md')).read_text(encoding='utf-8')
        reader_path = target / 'reader.json'
        if method == 'background':
            assert not list(unpacked.rglob('polyMesh')), 'background must not export a fluid case'
            assert not list(unpacked.rglob('*.solver.cm2d')), 'background must not masquerade as solver mesh'
            grid = next(unpacked.rglob('*.background.json'))
            assert next(unpacked.rglob('*.background.vtk')).is_file()
            reader = audit_background(grid)
            assert reader['full_domain_coverage'] and reader['classification_checked'] and reader['balanced']
            assert reader['counts'][1] > 0 and reader['counts'][2] > 0, 'solid and intersected cells retained'
            assert reader['solver_ready'] is False
            assert ui['background']['panelsHidden']
            assert set(ui['background']['rejected']) == {'runFlow', 'runThermal', 'runEuler'}
            reader.update(valid=True, cell_count=reader['cells'])
            reader_path.write_text(json.dumps(reader, indent=2), encoding='utf-8')
        else:
            case = next(unpacked.rglob('polyMesh')).parent.parent
            subprocess.run([sys.executable, str(ROOT / 'tools/verification/check_openfoam2d.py'),
                            str(case), '--report', str(reader_path)], check=True, stdout=subprocess.DEVNULL)
            reader = json.loads(reader_path.read_text(encoding='utf-8'))
        assert reader['valid'] and reader['cell_count'] == ui['layout']['previewCells']
        assert ui['layout']['exportedCells'] == (0 if method == 'background' else reader['cell_count'])
        if image:
            imported = json.loads(next(unpacked.rglob('image-import.json')).read_text(encoding='utf-8'))
            assert imported['physicalWidth'] == .2 and imported['outputUnits'] == 'm'
            assert imported['sourceSha256'] == hashlib.sha256(source.read_bytes()).hexdigest()
            assert next(unpacked.rglob('image-outline.png')).read_bytes()[:8] == b'\x89PNG\r\n\x1a\n'
        results.append({'case': name, 'cells': reader['cell_count'], 'independentReader': reader['valid'],
                        'layout': ui['layout'], 'sha256': hashlib.sha256(archive.read_bytes()).hexdigest()})
    (output / 'summary.json').write_text(json.dumps({'platform': sys.platform, 'cases': results,
        'externalCheckMesh': 'not_run', 'cfd': 'not_run'}, indent=2, ensure_ascii=False), encoding='utf-8')

if __name__ == '__main__':
    main()
