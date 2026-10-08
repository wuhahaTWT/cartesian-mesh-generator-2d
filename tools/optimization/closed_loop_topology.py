#!/usr/bin/env python3
"""Online sharp-wall feedback for native-2D density topology optimisation.

Porous adjoints propose material changes. Actual equal-area Cut-cell solves
measure their objective changes, and local discrepancy secants correct the next
material gradient. Topology-changing probes remain available independently of
the smooth local model. Only successfully solved sharp designs can be accepted.
No native adjoint, global optimum, or grid-independent accuracy is claimed.
"""
import argparse
from contextlib import redirect_stdout
from dataclasses import asdict, replace
import hashlib
import json
import math
import os
from pathlib import Path
import time

import numpy as np
from scipy.ndimage import gaussian_filter, label
from scipy.spatial import cKDTree
from scipy.optimize import minimize, LinearConstraint

from brinkman import Problem
from navier_stokes_brinkman import NavierStokesBrinkman
from optimize_flow import initial_design, metrics, snapshot, write_json
from compare_sharp_designs import prepare_design, flow_metrics
from topology_cases import problem_parameters
from engineering_baselines import distance_to_path
import native_flow as bridge

ROOT = Path(__file__).resolve().parents[2]


def topology_key(extraction):
    """Use extracted geometry, not density-cell labels, for model membership."""
    ports = extraction['portConnectivity']
    return json.dumps(dict(components=extraction['components'], holes=extraction['holes'],
        componentPorts=sorted(sorted(p) for p in ports['componentPorts'])), sort_keys=True)


def discrepancy_gradient(x, low, high, key, samples, design_mask):
    """Minimum-norm correction matching measured low/high objective secants.

    Each equation is s_i . correction = (H_i-H_0) - (L_i-L_0).
    Unit-normalized secants avoid favoring large perturbations. SVD discards
    near-dependent directions at sqrt(machine epsilon); vanishing material
    changes are also excluded before division, avoiding noise amplification at
    a porous stationary point. These are arithmetic safeguards, not CFD gates;
    this is a local regression, not an uncertainty estimator. Different sharp
    connectivity/holes are handled as explicit candidate branches, never mixed
    into this smooth correction. Failed native samples are excluded.
    """
    rows, rhs, ids = [], [], []
    for sample in samples:
        if sample.get('high') is None or sample.get('topology') != key:
            continue
        s = np.asarray(sample['x'])-x
        s = s.copy(); s[~design_mask] = 0
        length = float(np.linalg.norm(s))
        if np.max(np.abs(s)) <= math.sqrt(np.finfo(float).eps)*max(1., float(np.max(np.abs(x)))):
            continue
        rows.append(s/length)
        rhs.append(((sample['high']-high)-(sample['low']-low))/length)
        ids.append(sample['id'])
    correction = np.zeros_like(x)
    if not rows:
        return correction, dict(sampleIds=[], rank=0, correctionNorm=0., secantResidual=0.)
    matrix = np.asarray(rows)
    correction, _, rank, _ = np.linalg.lstsq(matrix, np.asarray(rhs),
                                           rcond=math.sqrt(np.finfo(float).eps))
    correction[~design_mask] = 0
    return correction, dict(sampleIds=ids, rank=int(rank),
        correctionNorm=float(np.linalg.norm(correction)),
        secantResidual=float(np.linalg.norm(matrix@correction-rhs)))


def material_probes(model, x, evaluation, beta, move, seed, round_index):
    """Adjoint, interface, and connectivity changes in the full material field.

    A shortest bridge between two density components is only a proposal. Its
    connectivity is measured again AFTER filtering, projection, and sharp-area
    matching. Passive collars and exact porous volume are preserved throughout.
    """
    proposals = [('porous-gradient', model.feasible_design(
        model.gradient_candidate(x, evaluation, beta, move), beta))]
    rho = evaluation.rho.reshape(model.ny, model.nx)
    yy, xx = np.mgrid[:model.ny, :model.nx]
    components, count = label(rho >= .5)
    if count > 1:
        # A thin nearest-component bridge can disappear in the equal-area
        # extraction. Also propose a whole shared corridor, generated from the
        # actual physical ports rather than prescribing a case-specific XY.
        # This is an explicit geometric search branch, not an adjoint result.
        def port_point(port):
            if port.side == 'left': return np.array([0., port.centre])
            if port.side == 'right': return np.array([model.problem.width, port.centre])
            if port.side == 'bottom': return np.array([port.centre, 0.])
            return np.array([port.centre, model.problem.height])
        inlet = [p for p in model.ports if p.role == 'inlet']
        outlet = [p for p in model.ports if p.role == 'outlet']
        def centre(ports):
            return np.average([port_point(p) for p in ports], axis=0,
                              weights=[p.width*p.peak for p in ports])
        left, right = centre(inlet), centre(outlet)
        fraction = .2+.1*(round_index % 3)
        a, b = left+fraction*(right-left), right-fraction*(right-left)
        paths = [np.array([port_point(p), a]) for p in inlet]+[np.array([b, port_point(p)]) for p in outlet]
        if np.linalg.norm(b-a) > 0:
            paths.append(np.array([a, b]))
        points = np.c_[(xx.ravel()+.5)*model.dx, (yy.ravel()+.5)*model.dy]
        # Opposing terminal barycentres can coincide; skip zero-length paths.
        paths = [p for p in paths if np.linalg.norm(p[1]-p[0]) > 0]
        signed = np.minimum.reduce([distance_to_path(points, p) for p in paths])
        lo, hi = 0., max(model.problem.width, model.problem.height)
        for _ in range(50):
            width = (lo+hi)/2
            routed = model.enforce_passive(np.clip(.5+(width-signed)/min(model.dx, model.dy), 0, 1))
            if model.volume(routed, beta)[0] > model.problem.volume_fraction:
                hi = width
            else:
                lo = width
        routed = model.enforce_passive(np.clip(.5+(lo-signed)/min(model.dx, model.dy), 0, 1))
        proposals.append(('merge-shared-corridor', model.feasible_design(routed, beta)))
        groups = sorted((np.argwhere(components == k) for k in range(1, count+1)),
                        key=len, reverse=True)
        a, b = groups[:2]
        scale = np.array([model.dy, model.dx])
        distances, indices = cKDTree(b*scale).query(a*scale)
        index = int(np.argmin(distances)); start, end = a[index]*scale, b[indices[index]]*scale
        points = np.stack((yy*model.dy, xx*model.dx), axis=-1)
        segment = end-start
        t = np.clip(np.sum((points-start)*segment, axis=-1)/max(float(segment@segment), 1e-30), 0, 1)
        distance = np.linalg.norm(points-start-t[..., None]*segment, axis=-1)
        # Proposal width is the existing physical filter length, not a new
        # certified minimum-channel-width or mesh-quality threshold.
        opening = np.exp(-.5*(distance/model.problem.filter_radius)**2).ravel()
        y = model.enforce_passive(np.clip(x+(1-x)*opening, 0, 1))
        proposals.append(('open-component-bridge', model.feasible_design(y, beta)))
    rng = np.random.default_rng(np.random.SeedSequence([seed, round_index]))
    noise = gaussian_filter(rng.standard_normal(rho.shape),
        sigma=(model.problem.filter_radius/model.dy, model.problem.filter_radius/model.dx))
    # Favor interface motion, while retaining material creation in solid areas.
    noise *= .25+4*rho*(1-rho)
    noise /= max(float(np.max(np.abs(noise))), np.finfo(float).tiny)
    for sign in (1, -1):
        y = model.enforce_passive(np.clip(x+sign*move*noise.ravel(), 0, 1))
        proposals.append(('interface-open' if sign == 1 else 'interface-close',
                          model.feasible_design(y, beta)))
    return proposals


def measured_basis(x, feedback, samples):
    selected = set(feedback['sampleIds'])
    directions = [s['x']-x for s in samples if s['id'] in selected]
    if not directions:
        return np.empty((len(x), 0))
    _, singular, vt = np.linalg.svd(np.asarray(directions), full_matrices=False)
    rank = int(np.sum(singular > math.sqrt(np.finfo(float).eps)*singular[0]))
    return vt[:rank].T


def corrected_material_search(model, x, evaluation, correction, q, beta, move, steps, basis=None):
    """Reoptimise L(x)+c.(x-x0) with the current measured discrepancy gradient.

    This inexpensive inner search creates a new full-field density design.
    The correction is held fixed during this local search; the native evaluator
    alone decides acceptance afterwards. The existing Armijo coefficient is
    reused from optimise_flow; it is not a new sharp-CFD acceptance threshold.
    """
    if basis is not None and basis.shape[1]:
        # Restrict exploitation to measured material directions. Unmeasured
        # directions are still explored through structural/interface probes.
        # This avoids treating a minimum-norm secant fit as a native gradient
        # over the entire high-dimensional material field.
        cache_z, cache_e = np.zeros(basis.shape[1]), evaluation
        def evaluate(z):
            nonlocal cache_z, cache_e
            if not np.array_equal(cache_z, z):
                y = model.enforce_passive(np.clip(x+basis@z, 0, 1))
                cache_e = model.evaluate(y, q, beta, 'total-pressure-power')
                cache_z = z.copy()
            return cache_e
        scale = max(abs(evaluation.objective), 1.)
        def objective(z):
            e = evaluate(z)
            return (e.objective+float(correction@(basis@z)))/scale
        def gradient(z):
            return basis.T@(evaluate(z).gradient+correction)/scale
        constraint = dict(type='ineq', fun=lambda z:model.problem.volume_fraction-evaluate(z).volume,
                          jac=lambda z:-basis.T@evaluate(z).volume_gradient)
        active = np.any(basis != 0, axis=1)
        lower, upper = np.maximum(-move, -x), np.minimum(move, 1-x)
        result = minimize(objective, np.zeros(basis.shape[1]), jac=gradient, method='SLSQP',
            constraints=[LinearConstraint(basis[active], lower[active], upper[active]), constraint],
            options=dict(maxiter=steps*10, ftol=1e-9))
        y = model.feasible_design(model.enforce_passive(np.clip(x+basis@result.x, 0, 1)), beta)
        return y, [dict(method='measured-material-subspace', dimensions=int(basis.shape[1]),
            optimizerStatus=int(result.status), optimizerMessage=str(result.message),
            optimizerIterations=int(result.nit), surrogateSearchConverged=bool(result.success),
            correctedMerit=float(result.fun*scale), materialChange=float(np.max(np.abs(y-x))))]
    y, e = x.copy(), evaluation
    history = []
    for _ in range(steps):
        updated = replace(e, gradient=e.gradient+correction)
        merit = e.objective+float(correction@(y-x))
        accepted = False
        step_move = move
        for _ in range(10):
            candidate = model.feasible_design(model.gradient_candidate(y, updated, beta, step_move), beta)
            predicted = float(updated.gradient@(candidate-y))
            if predicted >= 0:
                break
            try:
                trial = model.evaluate(candidate, q, beta, 'total-pressure-power')
                trial_merit = trial.objective+float(correction@(candidate-x))
                if trial_merit <= merit+1e-4*predicted:
                    accepted = True; break
            except (ArithmeticError, RuntimeError):
                pass
            step_move *= .5
        if not accepted:
            break
        y, e = candidate, trial
        history.append(dict(porousObjective=e.objective, correctedMerit=trial_merit,
                            materialChange=float(np.max(np.abs(y-x)))))
    return y, history


class NativeEvaluator:
    """One fixed equation/mesh contract for the entire feedback trajectory."""
    def __init__(self, root, model, q, beta, args):
        self.root, self.model, self.q, self.beta, self.args = root, model, q, beta, args
        self.samples, self.records, self.cache = [], [], {}
        self.started = time.monotonic()

    def evaluate(self, x, name, feedback=None, round_index=None):
        identity = hashlib.sha256(np.asarray(x, dtype='<f8').tobytes()).hexdigest()
        if identity in self.cache:
            return self.cache[identity]
        index = len(self.records)
        directory = self.root/'evaluations'/f'{index:03d}-{name}'
        directory.mkdir(parents=True)
        row = dict(id=index, name=name, round=round_index, directory=str(directory),
            materialSha256=identity, status='porous-analysis', high=None,
            feedback=feedback, elapsedSeconds=0.)
        self.records.append(row); self.save()
        started = time.monotonic()
        sample = dict(id=index, x=x.copy(), high=None, low=None, topology=None, row=row, e=None)
        self.cache[identity] = sample
        try:
            e = self.model.evaluate(x, self.q, self.beta, 'total-pressure-power')
            sample.update(e=e, low=e.objective)
            row['porous'] = metrics(self.model, x, e)
            snapshot(directory/'material.npz', self.model, x, e, self.q, self.beta)
            target = self.model.problem.width*self.model.problem.height*self.model.problem.volume_fraction
            sharp = prepare_design(directory/'sharp', directory/'material.npz', self.model, target, name)
            row['sharp'] = sharp
            sample['topology'] = row['topology'] = topology_key(sharp['extraction'])
            if not sharp['extraction']['portConnectivity']['allPortsCovered']:
                raise ValueError('sharp-area projection does not preserve all physical ports')
            row['status'] = 'native-solving'; self.save()
            a = self.args
            native_args = argparse.Namespace(directory=directory/'sharp', output=directory/'native',
                mesh_cli=a.mesh_cli, flow_cli=a.flow_cli, levels=[a.level], timeout=a.timeout,
                iterations=a.native_iterations, speed=a.speed, nu=1., tolerance=a.tolerance,
                small_alpha=a.small_alpha, padding_fraction=a.padding_fraction,
                mesh_only=False, component=None, reynolds=self.model.problem.reynolds,
                convergence='strict', velocity_relaxation=.2, steady_acceleration='anderson')
            with (directory/'native.log').open('w') as log, redirect_stdout(log):
                native = bridge.run(native_args)
            row['nativeSummary'] = str(directory/'native'/'summary.json')
            row['native'] = flow_metrics(native)
            if row['native'] is None or row['native']['dimensionlessTotalPower'] is None:
                row['status'] = 'native-rejected'
            else:
                high = row['native']['dimensionlessTotalPower']
                if not math.isfinite(high) or high <= 0:
                    raise ValueError('nonpositive or nonfinite sharp-wall objective')
                # Fixed prescribed port profiles make throughput a physical
                # constraint, not a degree of freedom the optimiser can exploit.
                expected = self.model.inflow*a.speed
                if abs(row['native']['inletFlux']-expected) > 1e-10*expected:
                    raise ValueError('native throughput differs from porous port contract')
                sample['high'] = row['high'] = high
                row['status'] = 'evaluated'
        except (ValueError, ArithmeticError, RuntimeError) as exc:
            row.update(status='evaluation-failed', issue=str(exc))
        row['elapsedSeconds'] = time.monotonic()-started
        self.samples.append(sample); self.save()
        print(f"eval={index} {name} status={row['status']} J_B={sample['low']} J_T={sample['high']}", flush=True)
        return sample

    def save(self):
        write_json(self.root/'evaluations.json', dict(evaluations=self.records,
            elapsedSeconds=time.monotonic()-self.started))


def evolve(model, x, q, beta, evaluator, rounds=4, move=.15, probes=2, seed=0,
           feedback_enabled=True, max_seconds=1800., material_steps=4, native_backtracks=1):
    """Online accepted-state loop; failed trials never replace the incumbent."""
    started = time.monotonic()
    incumbent = evaluator.evaluate(x, 'initial')
    if incumbent['high'] is None:
        raise RuntimeError('initial sharp-wall design failed; no accepted state available')
    report = dict(schema='cartmesh2d-sharp-feedback-v1', status='running',
        method='measured-subspace discrepancy optimisation + material/connectivity probes + native acceptance',
        feedbackEnabled=feedback_enabled, problem=asdict(model.problem),
        q=q, beta=beta, history=[], acceptedIds=[incumbent['id']],
        initialId=incumbent['id'], initialObjective=incumbent['high'],
        physicalAccuracyQualified=False, optimizationConverged=False,
        controls=dict(rounds=rounds, move=move, probes=probes, seed=seed,
                      maxSeconds=max_seconds, materialSteps=material_steps,
                      nativeBacktracks=native_backtracks), limits=[
            'Online new material designs; no exact native gradient or global-optimum claim.',
            'Objective improvements refer to one fixed native mesh/equation contract.',
            'Round/time limits are budgets, not stationarity or mesh-independence certificates.'])
    def save():
        report.update(finalId=incumbent['id'], finalObjective=incumbent['high'],
            relativeReduction=1-incumbent['high']/report['initialObjective'],
            elapsedSeconds=time.monotonic()-started,
            nativeEvaluations=len(evaluator.records),
            successfulNativeEvaluations=sum(r['high'] is not None for r in evaluator.records))
        write_json(evaluator.root/'summary.json', report)
    snapshot(evaluator.root/'accepted-design.npz', model, incumbent['x'], incumbent['e'], q, beta)
    save()
    for iteration in range(rounds):
        if time.monotonic()-started >= max_seconds:
            report['status'] = 'time-limit'; break
        x, e = incumbent['x'], incumbent['e']
        trials, step = [], dict(round=iteration, centreId=incumbent['id'],
            centreObjective=incumbent['high'], move=move, probeIds=[])
        report['history'].append(step); save()
        def evaluate_trial(proposal, name, feedback=None):
            trial = evaluator.evaluate(proposal, name, feedback, iteration)
            step.setdefault('trialIds', []).append(trial['id']); save()
            for retry in range(native_backtracks):
                if trial['high'] is not None or time.monotonic()-started >= max_seconds:
                    break
                step.setdefault('rejectedTrialIds', []).append(trial['id'])
                # Retry the MATERIAL perturbation, never change mesh quality,
                # physics or convergence controls to get an accepted score.
                proposal = model.feasible_design((x+proposal)/2, beta)
                trial = evaluator.evaluate(proposal, f'{name}-half-{retry+1}', feedback, iteration)
                step['trialIds'].append(trial['id']); save()
            return trial
        bank = material_probes(model, x, e, beta, move, seed, iteration)
        # Always include an interface probe: adjoint-only trials can be tiny at
        # a porous local optimum. Rotate the structural branch across rounds.
        selected = [bank[0]]
        structural = bank[1:]
        for j in range(max(0, probes-1)):
            selected.append(structural[(iteration+j) % len(structural)])
        for name, proposal in selected:
            if time.monotonic()-started >= max_seconds:
                break
            trial = evaluate_trial(proposal, name)
            trials.append(trial); step['probeIds'].append(trial['id']); save()
        if feedback_enabled and time.monotonic()-started < max_seconds:
            # Recent local measurements keep the secant model from expanding
            # indefinitely; branch filtering is done inside the regression.
            correction, feedback = discrepancy_gradient(x, e.objective, incumbent['high'],
                incumbent['topology'], evaluator.samples[-12:], model.design)
            feedback.update(centreId=incumbent['id'],
                lowGradientNorm=float(np.linalg.norm(e.gradient)),
                correctedGradientNorm=float(np.linalg.norm(e.gradient+correction)))
            step['feedback'] = feedback
            if feedback['rank']:
                basis = measured_basis(x, feedback, evaluator.samples[-12:])
                try:
                    proposal, inner = corrected_material_search(model, x, e, correction, q, beta,
                                                                move, material_steps, basis)
                    step['materialSearch'] = inner
                    step['correctedMaterialChange'] = float(np.max(np.abs(proposal-x)))
                    porous_proposal = bank[0][1]
                    step['changeFromPorousProposal'] = float(np.max(np.abs(proposal-porous_proposal)))
                    trial = evaluate_trial(proposal, 'feedback-gradient', feedback)
                    step['correctedId'] = trial['id']; trials.append(trial); save()
                except (ArithmeticError, RuntimeError, ValueError) as exc:
                    step['surrogateFailure'] = str(exc); save()
        successful = [t for t in trials if t['high'] is not None]
        best = min([incumbent]+successful, key=lambda t: t['high'])
        step.update(accepted=bool(best['high'] < incumbent['high']), chosenId=best['id'],
                    chosenObjective=best['high'])
        if step['accepted']:
            step['topologyChanged'] = best['topology'] != incumbent['topology']
            incumbent = best; report['acceptedIds'].append(best['id'])
            snapshot(evaluator.root/'accepted-design.npz', model, incumbent['x'], incumbent['e'], q, beta)
        else:
            move *= .5
        print(f"round={iteration} accepted={step['accepted']} chosen={best['id']} J_T={best['high']:.9g}", flush=True)
        save()
    else:
        report['status'] = 'round-limit'
    if report['status'] == 'running':
        report['status'] = 'time-limit'
    if time.monotonic()-started >= max_seconds:
        report['status'] = 'time-limit'
    report['feedbackDrivenAcceptedSteps'] = sum(
        h['accepted'] and h.get('correctedId') == h['chosenId'] and
        evaluator.records[h['chosenId']].get('feedback') is not None and
        h.get('changeFromPorousProposal', 0) > 0 for h in report['history'])
    report['acceptedTopologyChanges'] = sum(h.get('topologyChanged', False) for h in report['history'])
    save()
    snapshot(evaluator.root/'accepted-design.npz', model, incumbent['x'], incumbent['e'], q, beta)
    return report


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, help='existing porous optimisation directory; never modified')
    p.add_argument('--field', default='final.npz')
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--case', choices=['double-pipe', 'diffuser', 'elbow', 'four-terminal'], default='double-pipe')
    p.add_argument('--reynolds', type=float, default=0.)
    p.add_argument('--initialization', choices=['uniform', 'random', 'geometric', 'merged'], default='geometric')
    p.add_argument('--alpha-max', type=float, default=360000.)
    p.add_argument('--q', type=float, default=.1)
    p.add_argument('--beta', type=float, default=6.)
    p.add_argument('--rounds', type=int, default=4)
    p.add_argument('--move', type=float, default=.15)
    p.add_argument('--probes', type=int, default=2)
    p.add_argument('--material-steps', type=int, default=4,
                   help='cheap material optimisation steps per measured correction')
    p.add_argument('--native-backtracks', type=int, default=1,
                   help='halve a failed material perturbation; preserve the failed native attempt')
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--max-seconds', type=float, default=1800.)
    p.add_argument('--level', type=int, default=5)
    p.add_argument('--native-iterations', type=int, default=6000)
    p.add_argument('--timeout', type=float, default=120.)
    p.add_argument('--tolerance', type=float, default=1e-8)
    p.add_argument('--speed', type=float, default=.02)
    p.add_argument('--small-alpha', type=float, default=.25)
    p.add_argument('--padding-fraction', type=float, default=1/30)
    p.add_argument('--mesh-cli', type=Path, default=ROOT/'build/cartmesh2d_cli')
    p.add_argument('--flow-cli', type=Path, default=ROOT/'build/cartmesh2d_flow_cli')
    p.add_argument('--no-feedback', action='store_true', help='same probes/native acceptance, no gradient correction')
    return p


def run(args):
    if args.rounds < 1 or args.probes < 1 or args.material_steps < 1 or args.native_backtracks < 0:
        raise ValueError('rounds, probes and material steps must be positive; backtracks must be nonnegative')
    if not math.isfinite(args.move) or not 0 < args.move <= 1:
        raise ValueError('material move must be in (0, 1]')
    root = args.output.resolve()
    source_metadata = None
    if args.source:
        source = args.source.resolve(strict=True)
        original = json.loads((source/'summary.json').read_text())
        model = NavierStokesBrinkman(Problem(**original['problem']))
        field = (source/args.field).resolve(strict=True)
        with np.load(field) as data:
            if float(data['width']) != model.problem.width or float(data['height']) != model.problem.height:
                raise ValueError('source field dimensions differ from the recorded problem')
            x = model.enforce_passive(data['design'].ravel())
            q, beta = float(data['q']), float(data['beta'])
        source_metadata = dict(directory=str(source), field=str(field), sha256=bridge.sha(field))
    else:
        problem = Problem(**problem_parameters(args.case), case=args.case,
            filter_radius=.055, alpha_max=args.alpha_max, reynolds=args.reynolds)
        model = NavierStokesBrinkman(problem)
        q, beta = args.q, args.beta
        x = model.feasible_design(initial_design(model, args.initialization, args.seed), beta)
    args.mesh_cli = args.mesh_cli.resolve(strict=True)
    args.flow_cli = args.flow_cli.resolve(strict=True)
    root.mkdir(parents=True)
    os.environ.setdefault('MPLCONFIGDIR', str(root/'plot-cache'))
    write_json(root/'run.json', dict(source=source_metadata, problem=asdict(model.problem),
        q=q, beta=beta, nativeControls=dict(level=args.level, smallAlpha=args.small_alpha,
            paddingFraction=args.padding_fraction, iterations=args.native_iterations,
            speed=args.speed, tolerance=args.tolerance, timeout=args.timeout),
        sourceHashes={str(p.relative_to(ROOT)):bridge.sha(p) for p in
            (Path(__file__), Path(__file__).with_name('brinkman.py'),
             Path(__file__).with_name('navier_stokes_brinkman.py'),
             Path(__file__).with_name('compare_sharp_designs.py'), Path(bridge.__file__))},
        binaryHashes={str(p):bridge.sha(p) for p in (args.mesh_cli, args.flow_cli)}))
    evaluator = NativeEvaluator(root, model, q, beta, args)
    try:
        return evolve(model, x, q, beta, evaluator, args.rounds, args.move, args.probes,
            args.seed, not args.no_feedback, args.max_seconds, args.material_steps, args.native_backtracks)
    except RuntimeError as exc:
        if not (root/'summary.json').exists():
            write_json(root/'summary.json', dict(schema='cartmesh2d-sharp-feedback-v1',
                status='initial-evaluation-failed', issue=str(exc), acceptedIds=[],
                problem=asdict(model.problem), physicalAccuracyQualified=False))
        raise


if __name__ == '__main__':
    args = parser().parse_args()
    result = run(args)
    print(json.dumps({k: result[k] for k in ('status', 'initialObjective', 'finalObjective',
        'relativeReduction', 'feedbackDrivenAcceptedSteps', 'acceptedTopologyChanges')}, indent=2))
