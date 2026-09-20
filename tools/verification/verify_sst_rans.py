#!/usr/bin/env python3
"""Independent audit of the experimental coupled steady SST channel output.

This audit checks returned flow and turbulence fields, recomputes the SST
closure and both scalar transport balances, and delegates momentum-face
reconstruction to the existing native-flow auditor with independently derived
face viscosities.  It deliberately does not apply the constant-viscosity
Poiseuille benchmark.
"""
import argparse
import csv
import json
import math
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import verify_native_flow as native  # noqa: E402
from verify_sst_spatial import point_segment  # noqa: E402


def req(condition, message):
    if not condition:
        raise ValueError(message)


def num(value, label):
    try:
        x = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} is not numeric") from exc
    req(math.isfinite(x), f"{label} is not finite")
    return x


def rows(path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        req(reader.fieldnames is not None, f"missing CSV header: {path}")
        return list(reader), list(reader.fieldnames)


def close(a, b, absolute=1e-11, relative=1e-8):
    return abs(a - b) <= absolute + relative * max(abs(a), abs(b))


def coefficients(k, omega, nu, distance, strain, grad_k, grad_w):
    req(k >= 0 and omega > 0 and nu > 0 and distance > 0 and strain >= 0,
        "invalid SST point input")
    cross_base = 2 * .856 * (grad_k[0] * grad_w[0] + grad_k[1] * grad_w[1]) / omega
    cd = max(cross_base, 1e-10)
    turbulent = math.sqrt(k) / (.09 * omega * distance)
    viscous = 500 * nu / (distance * distance * omega)
    cross_bound = 4 * .856 * k / (cd * distance * distance)
    arg1 = min(max(turbulent, viscous), cross_bound)
    arg2 = max(2 * turbulent, viscous)
    f1 = math.tanh(min(arg1, 3.) ** 4)
    f2 = math.tanh(min(arg2, 5.) ** 2)
    sigma_k = f1 * .85 + (1 - f1)
    sigma_w = f1 * .5 + (1 - f1) * .856
    beta = f1 * .075 + (1 - f1) * .0828
    gamma = f1 * (5 / 9) + (1 - f1) * .44
    denominator = max(.31 * omega, strain * f2)
    req(denominator > 0, "invalid SST viscosity denominator")
    nut = .31 * k / denominator
    limited_strain2 = min(strain * strain, 10 * .09 * omega * denominator / .31)
    production_k = nut * limited_strain2
    production_w = gamma * limited_strain2
    cross = (1 - f1) * cross_base
    return {
        "F1": f1, "F2": f2, "nuT": nut,
        "Dk": nu + sigma_k * nut, "Dw": nu + sigma_w * nut,
        "sourceK": production_k, "sourceW": production_w + max(cross, 0.),
        "lossK": .09 * omega, "lossW": beta * omega + max(-cross, 0.) / omega,
    }


def read_artifacts(mesh_path, prefix):
    mesh = native.read_cm2d(mesh_path)
    measured = native.measure(mesh, 1e-10, 1e-9)
    req(not measured.issues, str(measured.issues))
    base = Path(prefix)
    meta = json.loads(Path(str(base) + ".json").read_text(encoding="utf-8"))
    crows, cfields = rows(Path(str(base) + ".cells.csv"))
    frows, ffields = rows(Path(str(base) + ".faces.csv"))
    hrows, hfields = rows(Path(str(base) + ".history.csv"))
    required_c = ("cell", "x", "y", "area", "u", "v", "p", "speed", "k", "omega",
                  "distance", "gradKx", "gradKy", "gradWx", "gradWy", "strain",
                  "F1", "F2", "nuT", "Dk", "Dw", "sourceK", "sourceW", "lossK", "lossW")
    required_f = ("face", "owner", "neighbour", "flux", "pressure", "advectionX", "advectionY",
                  "diffusionX", "diffusionY", "viscosity", "wall", "kBoundary", "omegaBoundary",
                  "kAdvection", "kDiffusion", "omegaAdvection", "omegaDiffusion")
    required_h = ("iteration", "momentumResidual", "continuity", "velocityChange", "pressureChange",
                  "kNorm", "omegaNorm", "kCellResidual", "omegaCellResidual", "turbulenceIterations")
    req(all(x in cfields for x in required_c), "incomplete SST-RANS cell schema")
    req(all(x in ffields for x in required_f), "incomplete SST-RANS face schema")
    req(all(x in hfields for x in required_h), "incomplete SST-RANS history schema")
    req(len(crows) == len(mesh.cells) and len(frows) == len(mesh.edges), "field row count mismatch")
    req(meta.get("case") in ("channel", "flatplate") and meta.get("model") == "SST-2003m" and
        meta.get("scope") == "coupled-steady-SST-2003m" and meta.get("converged") is True,
        "invalid SST-RANS metadata")
    req(meta.get("nu") == .001 and meta.get("speed") == 1 and meta.get("inletK") == .001 and
        meta.get("inletOmega") == 2 and meta.get("tolerance") == 1e-7,
        "changed SST-RANS physical configuration")
    if meta['case']=='flatplate':
        req(meta.get('flatPlateLeadingEdge')==.5 and
            meta.get('flatPlateTop') in ('pressure-farfield','symmetry'), 'changed flat plate probe configuration')
    else:
        req(meta.get('flatPlateLeadingEdge',0)==0 and
            meta.get('flatPlateTop','pressure-farfield')=='pressure-farfield', 'flat plate controls on channel')
    req(meta.get("cells") == len(mesh.cells) and meta.get("scalarRelativeTolerance") == 1e-9 and
        meta.get("scalarAbsoluteTolerance") == 1e-12 and meta.get("scalarCellTolerance") == 1e-9,
        "changed SST-RANS counts or scalar tolerances")
    req(meta.get("pressureConvention") == "p/rho (SST-2003m omits isotropic k stress)",
        "invalid pressure convention")
    req(meta.get("convection") == "upwind" and meta.get("viscousStress") == "symmetric" and
        meta.get("pressureDiscretization") == "shared-face-gauss",
        "unsupported flow discretization")
    req(len(hrows) == int(meta.get("iterations", -1)) and len(hrows) > 0, "history length mismatch")
    updates=meta.get('turbulenceUpdatesPerIteration',500) # legacy nested probe
    req(type(updates) is int and 1<=updates<=500,'invalid constitutive update limit')
    corrections=meta.get('scalarCorrectionsPerUpdate',0) # legacy complete frozen solve
    req(type(corrections) is int and corrections in (0,1),'unsupported scalar correction schedule')
    for expected, row in enumerate(hrows, 1):
        req(int(row["iteration"]) == expected, "history iteration ordering mismatch")
        for key in ("momentumResidual", "continuity", "velocityChange", "pressureChange",
                    "kNorm", "omegaNorm", "kCellResidual", "omegaCellResidual"):
            req(num(row[key], key) >= 0, f"history {key} is negative")
        req(0 <= int(row["turbulenceIterations"]) <= updates, "invalid turbulence iteration count")
    # Both files serialize the very same native value at 17 digits. Independent
    # reconstruction needs a roundoff allowance; duplicate export data does not.
    req(num(meta["momentumResidual"], "momentumResidual") ==
        num(hrows[-1]["momentumResidual"], "momentumResidual"),
        "summary momentum residual differs from final history")
    return mesh, measured, meta, crows, frows, hrows


def audit(mesh_path, prefix):
    mesh, m, meta, cells, faces, history = read_artifacts(Path(mesh_path), Path(prefix))
    geo = native.face_geometry(mesh, m)
    nu = .001
    case=meta['case']
    if case=='flatplate':
        req(all(close(a,b,1e-12,1e-10) for a,b in zip(m.bounds,(0.,0.,1.,1.))) and
            close(m.total_area,1.,1e-12,1e-10), 'flat plate diagnostic requires complete unit square')
    boundaries = native.flow_boundaries(mesh, m, case, 1., "reject",
                                        [num(r["flux"], "face flux") for r in faces],
                                        meta.get('flatPlateLeadingEdge'),meta.get('flatPlateTop','pressure-farfield'))
    def entering(fid):
        return boundaries['roles'][fid]=='inlet' or (
            boundaries['roles'][fid]=='farfield' and num(faces[fid]['flux'],'flux')<0)
    independent_continuity = native.continuity(
        mesh, m, [num(r["flux"], "face flux") for r in faces], 1., case, 1e-14, 1e-9)
    req(independent_continuity['nativeDefinitionContinuity'] < 1e-8 and
        independent_continuity['globalRelativeImbalance'] < 1e-8, 'independent continuity failed')
    req(close(num(meta['globalRelativeImbalance'], 'global continuity'),
              independent_continuity['globalRelativeImbalance'], 1e-14, 1e-6), 'global continuity summary mismatch')
    n = len(mesh.cells)
    k = [num(r["k"], "k") for r in cells]
    w = [num(r["omega"], "omega") for r in cells]
    u = [num(r["u"], "u") for r in cells]
    v = [num(r["v"], "v") for r in cells]
    req(all(x >= 0 for x in k) and all(x > 0 for x in w), "inadmissible turbulence field")
    kb = [0.] * len(mesh.edges); wb = [0.] * len(mesh.edges)
    kfixed = [False] * len(mesh.edges); wfixed = [False] * len(mesh.edges)
    walls = []
    for i, (edge, face) in enumerate(zip(mesh.edges, faces)):
        req(int(face["face"]) == i and int(face["owner"]) == edge.owner and
            int(face["neighbour"]) == edge.neighbour, "face ordering/topology mismatch")
        wall = int(face["wall"])
        expected_wall = int(edge.neighbour < 0 and boundaries["roles"][i] in ("wall", "lid"))
        req(wall == expected_wall, f"face {i}: wall role mismatch")
        if edge.neighbour < 0:
            q=num(face['flux'],'boundary flux')
            if wall: req(q==0, 'wall must be impermeable')
            elif boundaries['roles'][i]=='slip': req(q==0, 'symmetry boundary must be impermeable')
            elif boundaries['roles'][i]=='inlet':
                expected=sum(a*b for a,b in zip((boundaries['u'][i],boundaries['v'][i]),geo[i].area_vector))
                req(close(q,expected,1e-13,1e-11),'inlet flux mismatch')
            elif boundaries['roles'][i]=='outlet': req(q>=0, 'unconfigured outlet backflow')
        if wall:
            walls.append(i); kfixed[i] = wfixed[i] = True
            kb[i] = 0.
            normal_distance = sum((a-b)*s for a,b,s in zip(geo[i].centre, m.centroids[edge.owner], geo[i].area_vector)) / math.hypot(*geo[i].area_vector)
            req(normal_distance > 0, "invalid wall normal distance")
            wb[i] = 60 * nu / (.075 * normal_distance * normal_distance)
            req(close(num(face["kBoundary"], "wall k boundary"), kb[i], 1e-12, 1e-10) and
                close(num(face["omegaBoundary"], "wall omega boundary"), wb[i], 1e-10, 1e-10),
                "invalid resolved wall turbulence boundary")
        elif entering(i):
            kfixed[i] = wfixed[i] = True; kb[i] = .001; wb[i] = 2.
        req(close(num(face["kBoundary"], "k boundary"), kb[i], 1e-12, 1e-10) and
            close(num(face["omegaBoundary"], "omega boundary"), wb[i], 1e-12, 1e-10),
            f"face {i}: scalar boundary mismatch")
    req(walls, "channel has no resolved wall faces")
    kg = native.reconstruct_gradient(mesh, m, geo, k, kb, kfixed)
    wg = native.reconstruct_gradient(mesh, m, geo, w, wb, wfixed)
    # Flow gradients use the actual channel velocity constraints.
    ug = native.reconstruct_gradient(mesh, m, geo, u, boundaries["u"], boundaries["fixedU"])
    vg = native.reconstruct_gradient(mesh, m, geo, v, boundaries["v"], boundaries["fixedV"])
    strain = [math.hypot(math.sqrt(2.) * ug[i][0], math.sqrt(2.) * vg[i][1],
                         ug[i][1] + vg[i][0]) for i in range(n)]
    coeff = []; original_source=[]
    for i, row in enumerate(cells):
        req(int(row["cell"]) == i, "cell ordering mismatch")
        for key, value in (("x", m.centroids[i][0]), ("y", m.centroids[i][1]), ("area", m.areas[i]),
                           ('speed',math.hypot(u[i],v[i])), ("strain", strain[i])):
            actual = num(row[key], key)
            req(close(actual, value), f"cell {i}: {key} mismatch")
        for keys,expected in [(('gradKx','gradKy'),kg[i]),(('gradWx','gradWy'),wg[i])]:
            actual=tuple(num(row[key],key) for key in keys)
            roundoff=128*math.ulp(1.)*max(1.,math.hypot(*actual),math.hypot(*expected))
            req(all(close(a,b,1e-11+roundoff,1e-8) for a,b in zip(actual,expected)),f'cell {i}: gradient mismatch')
        d = num(row["distance"], "distance")
        req(d > 0, "invalid wall distance")
        expected_distance = min(
            point_segment(m.centroids[i], mesh.vertices[mesh.edges[f].v0],
                          mesh.vertices[mesh.edges[f].v1]) for f in walls)
        req(close(d, expected_distance, 1e-12, 1e-9), f"cell {i}: wall distance mismatch")
        c = coefficients(k[i], w[i], nu, expected_distance, strain[i], kg[i], wg[i]); coeff.append(c)
        # Original unsplit nonlinear source, evaluated at returned fields.
        beta=c['F1']*.075+(1-c['F1'])*.0828
        gamma=c['F1']*5/9+(1-c['F1'])*.44
        limited=min(strain[i]**2,10*.09*w[i]*max(.31*w[i],strain[i]*c['F2'])/.31)
        cross=(1-c['F1'])*2*.856/w[i]*sum(a*b for a,b in zip(kg[i],wg[i]))
        original_source.append((c['nuT']*limited-.09*w[i]*k[i],gamma*limited-beta*w[i]**2+cross))
        for key, value in c.items():
            req(close(num(row[key], key), value, 1e-11, 1e-8), f"cell {i}: stale SST closure {key}")

    # Independently derive the face viscosity and scalar diffusivities.
    face_nu = []
    for i, (edge, row) in enumerate(zip(mesh.edges, faces)):
        if i in walls:
            value = nu
        elif edge.neighbour >= 0:
            q = geo[i].neighbour_weight
            value = nu + (1-q) * coeff[edge.owner]["nuT"] + q * coeff[edge.neighbour]["nuT"]
        elif entering(i):
            ci = edge.owner
            inlet = coefficients(.001, 2., nu, num(cells[ci]["distance"], "distance"), strain[ci], kg[ci], wg[ci])
            value = nu + inlet["nuT"]
        else:
            value = nu + coeff[edge.owner]["nuT"]
        face_nu.append(value)
        req(close(num(row["viscosity"], "face viscosity"), value, 1e-12, 1e-8),
            f"face {i}: effective viscosity mismatch")

    # Check scalar constitutive face fluxes and steady cell balances.
    scalar_diagnostics = {}
    for name, values, grads, bc, fixed, diff_key, adv_key, source_key, loss_key, boundary_key in (
        ("k", k, kg, kb, kfixed, "Dk", "kAdvection", "sourceK", "lossK", "kBoundary"),
        ("omega", w, wg, wb, wfixed, "Dw", "omegaAdvection", "sourceW", "lossW", "omegaBoundary")):
        residual = [-m.areas[i]*original_source[i][0 if name=='k' else 1] for i in range(n)]
        diagonal = [m.areas[i] * coeff[i][loss_key] for i in range(n)]
        base = [m.areas[i] * coeff[i][source_key] for i in range(n)]
        norm_flux = 0.
        for fid, (edge, geom, row) in enumerate(zip(mesh.edges, geo, faces)):
            i = edge.owner; q = num(row["flux"], "flux")
            D = nu if fid in walls else (coeff[i][diff_key] if edge.neighbour < 0 else
                (1-geom.neighbour_weight)*coeff[i][diff_key] + geom.neighbour_weight*coeff[edge.neighbour][diff_key])
            if edge.neighbour >= 0:
                j = edge.neighbour; other = values[j] if q < 0 else values[i]
                adv = q * other
                gf = tuple((1-geom.neighbour_weight)*grads[i][a] + geom.neighbour_weight*grads[j][a] for a in (0, 1))
                diff = D * geom.transmissibility * (values[i] - values[j]) - D * sum(gf[a]*geom.correction[a] for a in (0, 1))
                diagonal[i] += D * geom.transmissibility + max(q, 0.)
                diagonal[j] += D * geom.transmissibility + max(-q, 0.)
            else:
                other = bc[fid] if q < 0 and fixed[fid] else values[i]
                adv = q * other
                if fixed[fid]:
                    diff = D * geom.transmissibility * (values[i] - bc[fid]) - D * sum(grads[i][a]*geom.correction[a] for a in (0, 1))
                    diagonal[i] += D * geom.transmissibility
                else:
                    diff = 0.
                    diagonal[i] += max(q, 0.)
            scalar_diffusion_key = "kDiffusion" if name == "k" else "omegaDiffusion"
            req(close(num(row[adv_key], adv_key), adv, 1e-10, 1e-8) and
                close(num(row[scalar_diffusion_key], "diffusion"), diff, 1e-10, 1e-8),
                f"face {fid}: {name} flux mismatch")
            flux = adv + diff; residual[i] += flux
            if edge.neighbour >= 0: residual[edge.neighbour] -= flux
            elif fixed[fid]:
                base[i] += D * geom.transmissibility * bc[fid]
            else:
                base[i] -= bc[fid] * math.hypot(*geom.area_vector)
            if edge.neighbour < 0 and q < 0:
                base[i] -= q * bc[fid]
            norm_flux = max(norm_flux, abs(flux))
        norm = math.sqrt(math.fsum(x*x for x in residual))
        max_cell = max(abs(x) / diagonal[i] for i, x in enumerate(residual))
        base_norm=math.sqrt(math.fsum(x*x for x in base))
        norm_target = float(meta["scalarAbsoluteTolerance"]) + float(meta["scalarRelativeTolerance"]) * base_norm
        # Independent geometry/arithmetic can differ at roundoff in sums of
        # opposing wall fluxes. Scale that budget with actual equation terms,
        # not a fixed allowance larger than the entire small-k stopping target.
        roundoff=128*math.ulp(1.)*max(base_norm,norm_flux*math.sqrt(n),1e-12)
        req(norm <= norm_target + roundoff and max_cell <= float(meta["scalarCellTolerance"]) + 1e-12,
            f"{name} independent residual exceeds scalar stopping gates")
        scalar_diagnostics[name] = {"residualNorm": norm, "normTarget": norm_target,
                                    "maxCellResidual": max_cell, "maxFaceFlux": norm_flux,
                                    "globalBalance": math.fsum(residual),"normRoundoffBudget":roundoff}
        history_norm = num(history[-1]["kNorm" if name == "k" else "omegaNorm"], "history norm")
        history_cell = num(history[-1]["kCellResidual" if name == "k" else "omegaCellResidual"], "history cell residual")
        req(close(history_norm, norm, roundoff, 1e-6) and close(history_cell, max_cell, 1e-12, 1e-6),
            f"{name} history does not describe independently reconstructed balance")

    # Reuse the independently implemented momentum auditor with a temporary
    # face-viscosity file derived above, never with native viscosity values.
    with tempfile.TemporaryDirectory(prefix="sst-rans-audit-") as temp:
        vf = Path(temp) / "face-viscosity.csv"
        with vf.open("w", newline="") as stream:
            out = csv.writer(stream); out.writerow(("face", "viscosity"))
            out.writerows((i, value) for i, value in enumerate(face_nu))
        flow_cells = [{"cell": int(r["cell"]), "x": num(r["x"], "x"), "y": num(r["y"], "y"),
                       "area": num(r["area"], "area"), "u": u[i], "v": v[i], "p": num(r["p"], "p"),
                       "speed": num(r["speed"], "speed")} for i, r in enumerate(cells)]
        flow_faces = [{"face": int(r["face"]), "owner": int(r["owner"]), "neighbour": int(r["neighbour"]),
                       "flux": num(r["flux"], "flux"), "pressure": num(r["pressure"], "pressure"),
                       "advectionX": num(r["advectionX"], "advectionX"), "advectionY": num(r["advectionY"], "advectionY"),
                       "diffusionX": num(r["diffusionX"], "diffusionX"), "diffusionY": num(r["diffusionY"], "diffusionY"),
                       "viscosity": num(r["viscosity"], "viscosity"), "wall": float(int(r["wall"]))} for r in faces]
        summary_fields = ("momentumResidual", "pressureForceX", "pressureForceY", "discreteForceX",
                          "discreteForceY", "reconstructedForceX", "reconstructedForceY", "forceX", "forceY",
                          "wallForceX", "wallForceY", "wallViscousForceX", "wallViscousForceY")
        req(all(field in meta for field in summary_fields), "SST-RANS JSON missing momentum summary fields")
        payload = {"case": case, "viscosityModel": "face-values", "viscosityFile": str(vf),
                   "flatPlateLeadingEdge":meta.get('flatPlateLeadingEdge'),
                   "flatPlateTop":meta.get('flatPlateTop','pressure-farfield'),
                   "convection": "upwind", "viscousStress": "symmetric", "outletBackflow": "reject",
                   "pressureBoundaryReconstruction": "one-sided-linear-2ring",
                   **{field: meta[field] for field in summary_fields}}
        momentum = native.reconstruct_momentum_audit(
            mesh, m, flow_cells, flow_faces, nu, 1., case, payload,
            pressure_boundary_reconstruction="one-sided-linear-2ring")
        req(all(value <= 5e-10 for value in momentum["maxFaceDeviation"].values()),
            "independent momentum face reconstruction differs")
        req(all(momentum["summaryDeviation"][field]["absolute"] <= 5e-10
                for field in summary_fields), "independent momentum summary differs")
        req(momentum["cellResidual"]["maxNormalized"] <= float(meta["tolerance"]),
            "independent momentum residual exceeds configured tolerance")
    req(num(history[-1]["momentumResidual"], "momentumResidual") <= float(meta["tolerance"]),
        "reported momentum residual exceeds configured tolerance")
    req(close(num(history[-1]["momentumResidual"], "momentumResidual"),
              momentum["cellResidual"]["maxNormalized"], 1e-12, 1e-6),
        "reported momentum residual differs from independent reconstruction")
    req(num(history[-1]["velocityChange"], "velocityChange") < float(meta["tolerance"]) and
        num(history[-1]["pressureChange"], "pressureChange") < float(meta["tolerance"]) and
        num(history[-1]["continuity"], "continuity") < 1e-8,
        "reported coupled convergence metrics exceed configured gates")
    req(close(num(history[-1]["continuity"], "continuity"),
              independent_continuity["nativeDefinitionContinuity"], 1e-14, 1e-6),
        "reported continuity differs from independent face balance")
    req(num(history[-1]["kCellResidual"], "kCellResidual") <= float(meta["scalarCellTolerance"]) and
        num(history[-1]["omegaCellResidual"], "omegaCellResidual") <= float(meta["scalarCellTolerance"]),
        "reported scalar cell residual exceeds configured tolerance")
    boundary_summary={}
    for role in ('inlet','outlet','wall','slip','farfield'):
        ids=[i for i,e in enumerate(mesh.edges) if e.neighbour<0 and boundaries['roles'][i]==role]
        fluxes=[num(faces[i]['flux'],'boundary flux') for i in ids]
        boundary_summary[role]={'faces':len(ids),'length':math.fsum(math.hypot(*geo[i].area_vector) for i in ids),
            'netOutwardFlux':math.fsum(fluxes),'inwardFlux':math.fsum(-q for q in fluxes if q<0),
            'outwardFlux':math.fsum(q for q in fluxes if q>0),'inflowFaces':sum(q<0 for q in fluxes)}
    wall_samples=[]
    if case=='flatplate':
        for i in walls:
            e=mesh.edges[i];length=math.hypot(*geo[i].area_vector)
            tau=num(faces[i]['diffusionX'],'wall tangential force')/length
            dn=abs(geo[i].centre[1]-m.centroids[e.owner][1])
            wall_samples.append({'face':i,'x':geo[i].centre[0],'xFromLeadingEdge':geo[i].centre[0]-.5,
                'length':length,'kinematicShear':tau,'Cf':2*tau,'yPlus':dn*math.sqrt(abs(tau))/nu})
        wall_samples.sort(key=lambda row:row['x'])
    return {"valid": True, "scope": meta["scope"], "case":case,
            "scalarCorrectionsPerUpdate":meta.get('scalarCorrectionsPerUpdate',0),
            "flatPlateTop":meta.get('flatPlateTop'),"boundarySummary":boundary_summary,
            "plateWallSamples":wall_samples,"wallSampleScope":"Current discrete wall traction and owner-centre y+; not an accuracy qualification",
            "cells": n, "iterations": meta["iterations"],
            "scalar": scalar_diagnostics, "momentum": momentum, "continuity": independent_continuity,
            "mesh": str(mesh_path), "prefix": str(prefix)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--prefix", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.mesh.resolve(), args.prefix.resolve())
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
