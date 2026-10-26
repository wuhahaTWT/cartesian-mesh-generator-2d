#!/usr/bin/env python3
"""Compile/run native SIMPLE source-choice ablations; no Python PDE solver.

Research only. Generated binaries and fields live under outputs/. The Gauss
pressure variant changes spatial discretization, not just iteration order.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[2]
GAUSS = r'''#pragma once
#include "FlowEquation2D.hpp"
namespace cartmesh2d::fv::solver_detail {
// Native 2-D Gauss-linear pressure choice. Unknown pressure traces use the
// owner value (zero-normal-gradient); prescribed pressure traces stay exact.
inline Vec referencePressureFaces(const FvMesh2D& m,const Vec& p,const Vec& bc,
                                 const std::vector<bool>& fixed,
                                 const detail::FlowGradientStencil2D& stencil) {
    (void)stencil;
    Vec values(m.faces.size());
    for(std::size_t id=0;id<values.size();++id) {
        const auto& f=m.faces[id];
        values[id]=f.neighbour?interpolate(f,p):fixed[id]?bc[id]:p[f.owner];
    }
    return values;
}
inline std::vector<Vector2D> referencePressureGradient(const FvMesh2D& m,
    const Vec& p,const Vec& bc,const std::vector<bool>& fixed,
    const detail::FlowGradientStencil2D& stencil) {
    return detail::conservativePressureGradient(m,referencePressureFaces(m,p,bc,fixed,stencil));
}
}
'''


def replace_once(text, before, after):
    if text.count(before) != 1:
        raise RuntimeError("Native source changed; review the ablation anchor: " + before[:90])
    return text.replace(before, after, 1)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', required=True, type=Path)
    parser.add_argument('--boundary', required=True, type=Path)
    parser.add_argument('--nu', required=True, type=float)
    parser.add_argument('--speed', required=True, type=float)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--convection', default='upwind', choices=['upwind', 'limited-linear', 'face-limited-linear'])
    parser.add_argument('--iterations', default=1500, type=int)
    parser.add_argument('--seconds', default=45, type=float)
    parser.add_argument('--tolerance', default=1e-8, type=float,
                        help='explicit research stopping tolerance; product defaults are unchanged')
    parser.add_argument('--build-only', action='store_true', help='compile variants without running CFD')
    parser.add_argument('--coupling', default='simple', choices=['simple', 'simplec'])
    parser.add_argument('--steady-acceleration', default='none', choices=['none', 'anderson'])
    parser.add_argument('--anderson-history', default=4, type=int)
    parser.add_argument('--anderson-backtracking', default='on', choices=['on', 'off'],
                        help='research-only source ablation of the two fallback halvings')
    parser.add_argument('--linear-policy', default='strict', choices=['strict', 'adaptive'])
    parser.add_argument('--pressure-trace', default='gauss-owner',
                        choices=['gauss-owner', 'linear-extrapolated', 'skew-owner', 'linear-wall'],
                        help='Gauss gradient with unskewed/owner, original LS face traces, '
                             'LS internal skew only, or LS wall extrapolation only')
    parser.add_argument('--variants', nargs='+', choices=['baseline', 'velocity', 'pressure', 'both'],
                        default=['baseline', 'velocity', 'pressure', 'both'])
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    mesh, boundary = args.mesh.resolve(), args.boundary.resolve()
    build = ROOT / 'build'
    windows = sys.platform == 'win32'
    compiler = shutil.which('cl') if windows else '/usr/bin/clang++' if sys.platform == 'darwin' else '/usr/bin/c++'
    if not compiler:
        raise RuntimeError('Run from a Visual Studio x64 developer command prompt (cl.exe is required)')
    app_object = build / 'cartmesh2d_flow_cli.dir/Release/cartmesh2d_flow_cli.obj' if windows else build / 'CMakeFiles/cartmesh2d_flow_cli.dir/apps/cartmesh2d_flow_cli.cpp.o'
    libraries = [build / 'Release/cartmesh2d_fv.lib', build / 'Release/cartmesh2d.lib'] if windows else [build / 'libcartmesh2d_fv.a', build / 'libcartmesh2d.a']
    product_binary = build / 'Release/cartmesh2d_flow_cli.exe' if windows else build / 'cartmesh2d_flow_cli'
    for path in [mesh, boundary, app_object, *libraries, product_binary]:
        if not path.is_file():
            raise FileNotFoundError(path)
    sources = ['src/fv/Incompressible2D.cpp', 'src/fv/FlowEquation2D.cpp',
               'src/fv/FlowPressure2D.hpp', 'src/fv/FlowFlux2D.cpp', 'src/fv/FlowForces2D.cpp',
               'src/fv/FlowCoupled2D.cpp', 'src/fv/FlowSolverDetail2D.hpp',
               'include/cartmesh2d/fv/Incompressible2D.hpp',
               'include/cartmesh2d/fv/detail/FlowFaceOperators2D.hpp',
               'apps/cartmesh2d_flow_cli.cpp', 'tools/benchmarks/simple_source_compare.py']
    manifest = {
        'format': 'cartmesh2d-simple-source-comparison-v1',
        'sourceHead': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'sourceSha256': {p: digest(ROOT / p) for p in sources},
        'linkedLibrarySha256': {p.name: digest(p) for p in libraries},
        'inputSha256': {str(p): digest(p) for p in [mesh, boundary]},
        'compiler': (subprocess.run([compiler], capture_output=True, text=True).stderr.splitlines()[0]
                     if windows else subprocess.check_output([compiler, '--version'], text=True).splitlines()[0]),
        'scope': 'native 2-D SIMPLE ablations; external OpenFOAM/MOOSE solvers are NOT executed',
        'controls': {'velocityRelaxation': .6, 'pressureRelaxation': .25,
                     'pressureCorrectionPasses': 4, 'tolerance': args.tolerance,
                     'maxIterations': args.iterations, 'seconds': args.seconds,
                     'convection': args.convection, 'nu': args.nu, 'speed': args.speed},
        'buildCommands': [], 'runs': []}
    manifest['controls'].update(coupling=args.coupling, steadyAcceleration=args.steady_acceleration,
                               andersonHistory=args.anderson_history, linearPolicy=args.linear_policy,
                               pressureTrace=args.pressure_trace, andersonBacktracking=args.anderson_backtracking)
    record = out / 'comparison.json'

    def save():
        record.write_text(json.dumps(manifest, indent=2) + '\n')

    def build_command(command):
        command = list(map(str, command))
        manifest['buildCommands'].append(command)
        save()
        with (out / 'build.log').open('a') as log:
            subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)

    def compile_source(path):
        obj = path.with_suffix('.obj' if windows else '.o')
        if windows:
            build_command([compiler, '/nologo', '/std:c++20', '/O2', '/DNDEBUG', '/MT', '/EHsc',
                           '/utf-8', '/DNOMINMAX', '/DWIN32_LEAN_AND_MEAN',
                           '/I' + str(ROOT / 'include'), '/I' + str(ROOT / 'src/fv'),
                           '/I' + str(ROOT / 'apps'), '/I' + str(out), '/c', path, '/Fo' + str(obj)])
        else:
            build_command([compiler, '-std=c++20', '-O3', '-DNDEBUG', '-I' + str(ROOT / 'include'),
                           '-I' + str(ROOT / 'src/fv'), '-I' + str(ROOT / 'apps'),
                           '-I' + str(out), '-c', path, '-o', obj])
        return obj

    pressure_choice = GAUSS
    if args.pressure_trace != 'gauss-owner':
        start = pressure_choice.index('    (void)stencil;')
        end = pressure_choice.index('\n    return values;', start) + len('\n    return values;')
        trace = '    auto values=detail::pressureFaceValues(m,p,stencil.apply(p,bc),bc,fixed);\n'
        if args.pressure_trace == 'skew-owner':
            trace += ('    for(std::size_t id=0;id<values.size();++id)\n'
                      '        if(!m.faces[id].neighbour && !fixed[id])values[id]=p[m.faces[id].owner];\n')
        elif args.pressure_trace == 'linear-wall':
            trace += ('    for(std::size_t id=0;id<values.size();++id)\n'
                      '        if(m.faces[id].neighbour)values[id]=interpolate(m.faces[id],p);\n')
        pressure_choice = (pressure_choice[:start] +
            trace + '    return values;' + pressure_choice[end:])
        pressure_choice = pressure_choice.replace('// Native 2-D Gauss-linear pressure choice. Unknown pressure traces use the\n'
            '// owner value (zero-normal-gradient); prescribed pressure traces stay exact.',
            '// Research Gauss gradient with ' + args.pressure_trace + ' pressure traces.')
    (out / 'reference-pressure-choice.hpp').write_text(pressure_choice)
    equation = (ROOT / 'src/fv/FlowEquation2D.cpp').read_text()
    equation = replace_once(equation, '#include "FlowEquation2D.hpp"', '#include "reference-pressure-choice.hpp"')
    equation = replace_once(equation,
        '    const auto gradient = pressure_.apply(pressure, boundaryPressure);\n'
        '    return detail::conservativePressureGradient(mesh_,\n'
        '        detail::pressureFaceValues(mesh_, pressure, gradient, boundaryPressure, fixedPressure));',
        '    return referencePressureGradient(mesh_,pressure,boundaryPressure,fixedPressure,pressure_);')
    equation = replace_once(equation, 'result.gp = pressure_.apply(field.p, boundaryPressure);',
                            'result.gp = referencePressureGradient(mesh_,field.p,boundaryPressure,boundary.fixedP,pressure_);')
    equation = replace_once(equation,
        'result.pressureFaces = detail::pressureFaceValues(mesh_, field.p, result.gp,\n'
        '        boundaryPressure, boundary.fixedP);',
        'result.pressureFaces = referencePressureFaces(mesh_,field.p,boundaryPressure,boundary.fixedP,pressure_);')
    equation_path = out / 'gauss-pressure-equation.cpp'
    equation_path.write_text(equation)
    pressure = (ROOT / 'src/fv/FlowPressure2D.hpp').read_text()
    pressure = replace_once(pressure, '#include "FlowEquation2D.hpp"', '#include "reference-pressure-choice.hpp"')
    pressure = replace_once(pressure, 'pressureGradientStencil.apply(pc,zeros)',
                            'referencePressureGradient(m,pc,zeros,b.fixedP,pressureGradientStencil)')
    (out / 'gauss-pressure-correction.hpp').write_text(pressure)
    # Avoid inheriting the product's one-sided pressure-reconstruction label.
    app_path = out / 'research-flow-cli.cpp'
    app_path.write_text(replace_once((ROOT / 'apps/cartmesh2d_flow_cli.cpp').read_text(),
        'one-sided-linear-adaptive', 'research-gauss-' + args.pressure_trace))
    pressure_objects = None
    for variant in args.variants:
        binary = product_binary
        if variant != 'baseline' or args.anderson_backtracking == 'off':
            gauss = variant in ['pressure', 'both']
            source = (ROOT / 'src/fv/Incompressible2D.cpp').read_text()
            if args.anderson_backtracking == 'off':
                source = replace_once(source, 'backtrack<3;', 'backtrack<1;')
            if variant in ['velocity', 'both']:
                for component in ['u', 'v']:
                    source = replace_once(source,
                        f'r.{component}[i]=finite(r.{component}[i]-pressureResponse[i]*gc[i].{ "x" if component == "u" else "y"});',
                        f'r.{component}[i]=finite(r.{component}[i]-c.pressureRelaxation*pressureResponse[i]*gc[i].{ "x" if component == "u" else "y"});')
            if gauss:
                source = replace_once(source, '#include "FlowPressure2D.hpp"', '#include "gauss-pressure-correction.hpp"')
            path = out / (variant + '-incompressible.cpp')
            path.write_text(source)
            objects = [compile_source(path)]
            if gauss and pressure_objects is None:
                pressure_objects = [compile_source(equation_path), compile_source(app_path)]
            selected_app = pressure_objects[1] if gauss else app_object
            if gauss:
                objects.append(pressure_objects[0])
            binary = out / ('flow-' + variant + ('.exe' if windows else ''))
            if windows:
                build_command([compiler, '/nologo', selected_app, *objects, *libraries, '/Fe' + str(binary),
                               '/link', '/MANIFEST:EMBED', '/MANIFESTINPUT:' + str(ROOT / 'apps/cartmesh2d.exe.manifest')])
            else:
                build_command([compiler, selected_app, *objects, *libraries,
                               *(['-framework', 'Accelerate'] if sys.platform == 'darwin' else []), '-o', binary])
        if args.build_only:
            manifest.setdefault('builtVariants', {})[variant] = {
                'binary': str(binary), 'sha256': digest(binary)}
            save()
            print('built', variant, flush=True)
            continue
        prefix = out / variant
        command = [str(binary), '--mesh', str(mesh), '--boundary', str(boundary), '--case', 'custom',
                   '--nu', str(args.nu), '--speed', str(args.speed), '--convection', args.convection,
                   '--coupling', args.coupling, '--velocity-relaxation', '.6', '--pressure-relaxation', '.25',
                   '--pressure-corrections', '4', '--linear-policy', args.linear_policy, '--tolerance', str(args.tolerance),
                   '--steady-acceleration', args.steady_acceleration, '--anderson-history', str(args.anderson_history),
                   '--max-iterations', str(args.iterations), '--profile', '--output', str(prefix)]
        run = {'variant': variant, 'command': command, 'binarySha256AtLaunch': digest(binary),
               'changesSpatialDiscretization': variant in ['pressure', 'both']}
        manifest['runs'].append(run)
        save()
        start = time.monotonic()
        with prefix.with_suffix('.log').open('w') as log:
            try:
                result = subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, timeout=args.seconds)
                run['exitCode'] = result.returncode
            except subprocess.TimeoutExpired:
                run['exitCode'] = 'timeout'
        run['wallSeconds'] = time.monotonic() - start
        summary = prefix.with_suffix('.json')
        if summary.exists():
            run['summary'] = json.loads(summary.read_text())
        run['lastLogLine'] = prefix.with_suffix('.log').read_text().splitlines()[-1]
        save()
        print(variant, run['exitCode'], run.get('summary', {}).get('iterations'), flush=True)


if __name__ == '__main__':
    main()
