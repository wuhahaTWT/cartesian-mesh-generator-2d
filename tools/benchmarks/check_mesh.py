"""Run the actual OpenFOAM checkMesh utility, without project quality thresholds."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import time

def mesh_hashes(case):
    return {name: hashlib.sha256((case/'constant/polyMesh'/name).read_bytes()).hexdigest()
            for name in ('points', 'faces', 'owner', 'neighbour', 'boundary')}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case', type=Path)
    parser.add_argument('--docker', help='Docker executable; used when checkMesh is not on PATH')
    parser.add_argument('--image', default='opencfd/openfoam-run:2606')
    args = parser.parse_args()
    case = args.case.resolve()
    hashes = mesh_hashes(case)
    local = shutil.which('checkMesh')
    docker = args.docker or shutil.which('docker')
    if not local and not docker and os.name == 'nt':
        candidate = Path('C:/Program Files/Docker/Docker/resources/bin/docker.exe')
        if candidate.is_file():
            docker = str(candidate)
    if not local and not docker:
        raise SystemExit('OpenFOAM checkMesh/Docker not found; no project check substituted.')
    results = []
    for name, options, log_name in [
        ('standard', [], 'checkMesh.log'),
        ('extended', ['-allGeometry', '-allTopology'], 'checkMesh-extended.log')
    ]:
        if local:
            command = [local, '-case', str(case), *options]
        else:
            # The script is constant; the case path is only a volume argument.
            script = ('source /usr/lib/openfoam/openfoam2606/etc/bashrc >/dev/null 2>&1'
                      ' && exec checkMesh -case /case "$@"')
            command = [docker, 'run', '--rm', '--network', 'none',
                       '-v', str(case)+':/case', '--entrypoint', '/bin/bash', args.image,
                       '-lc', script, 'checkMesh', *options]
        started = time.monotonic()
        with (case/log_name).open('w', encoding='utf-8') as log:
            completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
        output = (case/log_name).read_text(encoding='utf-8', errors='replace')
        passed = completed.returncode == 0 and 'Mesh OK.' in output
        outcome = 'passed' if passed else (
            'reported_issues' if re.search(r'Failed\s+\d+\s+mesh checks', output)
            else 'tool_error')
        version = re.search(r'Version:\s*([^|\r\n]+)', output)
        results.append(dict(name=name, outcome=outcome, exit_code=completed.returncode,
                            seconds=time.monotonic()-started, command=command,
                            version=version.group(1).strip() if version else None,
                            log=log_name))
        print(name, outcome, 'exit='+str(completed.returncode), flush=True)
    report = dict(software='OpenFOAM checkMesh', criteria='tool defaults; no project meshQualityDict',
                  case=str(case), mesh_sha256=hashes,
                  mesh_unchanged=hashes == mesh_hashes(case), checks=results)
    (case/'checkMesh.json').write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n',
                                      encoding='utf-8')
    return 0 if report['mesh_unchanged'] and all(c['outcome']=='passed' for c in results) else 1

if __name__ == '__main__':
    raise SystemExit(main())
