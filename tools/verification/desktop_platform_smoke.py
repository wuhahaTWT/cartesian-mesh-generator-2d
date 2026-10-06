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

ROOT = Path(__file__).resolve().parents[2]

def compressible_controls(args, output, workspace, env):
    """Check actual packaged processes and saved bytes, without another PDE audit."""
    results = []
    summary = {'scope': 'packaged compressible control/reopen protocol; no physical qualification or native file-dialog coverage',
               'platform': sys.platform, 'cases': results, 'status': 'running'}
    def save():
        (output/'summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding='utf-8')
    def digest(file):
        return hashlib.sha256(file.read_bytes()).hexdigest()
    def launch(name, extra):
        target = workspace/name; target.mkdir()
        screenshot = target/'app.png'; archive = target/'result.zip'
        command = [str(args.app.resolve()), *args.electron_arg, '--smoke=rectangle', '--method=cutcell', '--control=manual',
                   '--wall-relative-size=.25', '--background-relative-size=.25', '--reference-length=1',
                   '--band-cells=1', '--padding-relative-size=.25', '--out='+str(target/'cases'),
                   '--shot='+str(screenshot), '--export='+str(archive), *extra]
        (target/'command.json').write_text(json.dumps(command, indent=2, ensure_ascii=False), encoding='utf-8')
        with (target/'app.log').open('w', encoding='utf-8') as log:
            completed = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=240)
        if completed.returncode:
            raise RuntimeError(f'{name}: App failed; see {target/"app.log"}')
        ui = json.loads(Path(str(screenshot)+'.json').read_text(encoding='utf-8'))
        assert ui['bundledChineseFontLoaded'] and ui['euler'], 'packaged renderer result missing'
        unpacked = target/'unpacked'
        with zipfile.ZipFile(archive) as package:
            assert package.testzip() is None
            package.extractall(unpacked)
        meshes = list(unpacked.rglob('*.solver.cm2d'))
        assert len(meshes) == 1, 'expected one final rectangle solver mesh'
        mesh = meshes[0]; euler = ui['euler']
        assert euler['fieldCells'] == euler['summary']['cells'] == ui['layout']['previewCells']
        assert next(unpacked.rglob('euler-preview.png')).read_bytes().startswith(b'\x89PNG\r\n\x1a\n')
        manifest = mesh.parent/euler['manifest']; record = json.loads(manifest.read_text(encoding='utf-8'))
        checkpoint = manifest.parent/'euler.checkpoint'
        assert record['status'] == 'complete' and record['acceptedTime'] == euler['summary']['time']
        assert digest(checkpoint) == record['checkpointSha256']
        results.append({'case': name, 'cells': euler['summary']['cells'], 'checks': euler['checks'],
                        'meshSha256': digest(mesh), 'checkpointSha256': digest(checkpoint), 'archiveSha256': digest(archive),
                        'acceptedTime': record['acceptedTime']})
        save()
        return ui, mesh, checkpoint
    save()
    try:
        ui, mesh, _ = launch('controls', ['--euler-controls=true'])
        checks = ui['euler']['checks']
        for key in ['implicitReservoirReachedNative', 'physicalRestartLocked', 'failedBudgetPreservedComplete',
                    'implicitResumeChecked', 'integratorChangeClearedStaleResult', 'stageGuardedReachedNative',
                    'initialPerturbationReachedNative', 'steadyThreeGatesChecked', 'earlySteadyProgress']:
            assert checks[key] is True, key
        seed_manifest = mesh.parent/checks['reopenSeedManifest']
        seed = json.loads(seed_manifest.read_text(encoding='utf-8'))
        seed_checkpoint = seed_manifest.parent/'euler.checkpoint'
        assert seed['status'] == 'complete' and seed['acceptedTime'] == checks['reopenSeedTime']
        assert digest(seed_checkpoint) == seed['checkpointSha256']
        fixture = workspace/'unfinished-state'; fixture.mkdir()
        shutil.copyfile(seed_checkpoint, fixture/'euler.checkpoint')
        unfinished = dict(seed, status='running')
        del unfinished['checkpointSha256']; del unfinished['acceptedTime']
        fixture_manifest = fixture/'desktop-state.json'
        fixture_manifest.write_text(json.dumps(unfinished, indent=2), encoding='utf-8')
        original_manifest_sha = digest(fixture_manifest); original_state_sha = digest(fixture/'euler.checkpoint')
        summary['unfinishedFixture'] = 'Copied real completed native field with an unfinalized manifest; this is not a simulated whole-App crash.'
        target_time = 2*seed['acceptedTime']
        options = ['--euler-reopen-manifest='+str(fixture_manifest), '--euler-reopen-end='+str(target_time)]
        reopened = []
        for name in ['reopen', 'reopen-repeat']:
            restored, regenerated_mesh, checkpoint = launch(name, options)
            check = restored['euler']['checks']
            for key in ['reopenedManifest', 'unfinishedLabelChecked', 'physicalControlsRestoredAndLocked', 'resumedTargetReached']:
                assert check[key] is True, key
            assert check['sourceTime'] == seed['acceptedTime'] and check['sourceCheckpointSha256'] == original_state_sha
            assert mesh.read_bytes() == regenerated_mesh.read_bytes(), 'regenerated mesh differs in fresh App process'
            assert restored['euler']['summary']['time'] == target_time
            reopened.append(checkpoint)
        assert reopened[0].read_bytes() == reopened[1].read_bytes(), 'independent fresh-process resumes differ'
        assert digest(fixture_manifest) == original_manifest_sha and digest(fixture/'euler.checkpoint') == original_state_sha
        summary.update(status='passed', independentResumeCheckpointsIdentical=True, sourceFilesUnchanged=True)
        save()
    except Exception as error:
        summary.update(status='failed', error=str(error)); save(); raise

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--app', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--electron-arg', action='append', default=[])
    parser.add_argument('--scope', choices=['full', 'compressible'], default='full')
    args = parser.parse_args()
    output = args.out.resolve()
    output.mkdir(parents=True, exist_ok=args.scope != 'compressible')
    results = []
    # Exercise UTF-8 plus spaces through Node, native argv/filesystem and ZIP.
    workspace = output / '中文 路径'
    workspace.mkdir(exist_ok=True)
    temporary = workspace / 'temp'; temporary.mkdir(exist_ok=True)
    env = dict(os.environ, TEMP=str(temporary), TMP=str(temporary), TMPDIR=str(temporary))
    if args.scope == 'compressible':
        compressible_controls(args, output, workspace, env)
        return
    # Retained legacy full-scope readers are not invoked by the compressible
    # protocol lane, whose field validation stays in the native/App pathway.
    from check_background_grid import audit as audit_background
    from verify_euler import audit as audit_euler
    cases = [('png', 'raster-input-L.png', 'cutcell', 'circle'),
             ('jpg', 'raster-input-L.jpg', 'cutcell', 'circle'),
             ('hybrid', None, 'hybrid', 'circle'),
             ('background-uniform', None, 'background', 'circle'),
             ('background-adaptive', None, 'background', 'circle'),
             ('euler-sealed', None, 'cutcell', 'rectangle'),
             ('euler-external', None, 'cutcell', 'circle')]
    for name, image, method, sample in cases:
        target = workspace / name; target.mkdir(exist_ok=True)
        screenshot = target / 'app.png'
        archive = target / 'result.zip'
        command = [str(args.app.resolve()), *args.electron_arg, '--smoke=' + sample,
                   '--method=' + method, '--out=' + str(target / 'cases'),
                   '--export=' + str(archive), '--shot=' + str(screenshot)]
        if name.startswith('euler-'):
            # Short development cases: these properties are test controls, not an air model.
            # Match the bounded rectangle/circle meshes from the branch App evidence.
            command += ['--control=manual', '--band-cells=3',
                        '--wall-relative-size=' + ('.0625' if sample == 'rectangle' else '.015625'),
                        '--background-relative-size=' + ('.09375' if sample == 'rectangle' else '21'),
                        '--padding-relative-size=' + ('.25' if sample == 'rectangle' else '10'),
                        '--small-alpha=' + ('.1' if sample == 'rectangle' else '.15')]
            command += ['--euler=true', '--euler-case=' + name.removeprefix('euler-'),
                        '--euler-flux=hllc', '--euler-order=2', '--euler-wall-gradient=quadratic',
                        '--euler-wall-model=no-slip', '--euler-conductivity=100',
                        '--euler-wall-thermal=temperature', '--euler-wall-value=400', '--euler-u=50',
                        '--euler-viscosity=' + ('2' if sample == 'rectangle' else '.02'),
                        '--euler-end-time=' + ('.001' if sample == 'rectangle' else '.00005')]
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
        euler_audits = []
        if name.startswith('euler-'):
            euler = ui['euler']
            checks = euler['checks']
            for check in ['wallGradientControlsReachedNative', 'wallGradientChangeClearedStaleResult',
                          'viscousControlsReachedNative', 'thermalControlsReachedNative',
                          'repeatedFieldsAndHistoryIdentical', 'resumeChecked',
                          'failedBudgetPreservedComplete', 'cancelResumeChecked', 'allFiveFieldMaps']:
                assert checks[check] is True, check
            assert euler['summary']['wallGradient'] == 'quadratic'
            assert euler['summary']['dynamicViscosity'] > 0 and euler['summary']['thermalConductivity'] > 0
            # Audit the exported complete result using only exported geometry and fields.
            mesh = next(unpacked.rglob('*.solver.cm2d'))
            checkpoint = mesh.parent / euler['files']['.checkpoint']
            assert checkpoint.is_file(), 'export lost the selected compressible checkpoint'
            audit = audit_euler(mesh, checkpoint.with_suffix(''))
            assert audit['valid'] is True
            euler_audits.append(audit)
            (target / 'euler-audit.json').write_text(json.dumps(audit, indent=2), encoding='utf-8')
        results.append({'case': name, 'cells': reader['cell_count'], 'independentReader': reader['valid'],
                        'eulerAudits': euler_audits, 'eulerChecks': ui.get('euler', {}).get('checks') if ui.get('euler') else None,
                        'layout': ui['layout'], 'sha256': hashlib.sha256(archive.read_bytes()).hexdigest()})
    (output / 'summary.json').write_text(json.dumps({'platform': sys.platform, 'cases': results,
        'externalCheckMesh': 'not_run', 'cfd': 'two short compressible wall cases; not physical qualification'}, indent=2, ensure_ascii=False), encoding='utf-8')

if __name__ == '__main__':
    main()
