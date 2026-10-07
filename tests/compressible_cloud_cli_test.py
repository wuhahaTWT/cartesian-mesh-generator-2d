"""Native generated-grid regression of termination, reservoir and restart semantics."""
import json
from pathlib import Path
import signal
import subprocess
import sys
import tempfile

mesh_cli, flow_cli = sys.argv[1:3]
with tempfile.TemporaryDirectory(prefix="cartmesh-cloud-") as folder:
    root = Path(folder)
    xy = root/"channel.xy"
    xy.write_text("0 0\n.0004 0\n.0004 .0001\n0 .0001\n")
    subprocess.run([mesh_cli, str(xy), str(root/"mesh"), "5", ".125", ".1", "interior", str(root/"foam"), "4", "0"], check=True, capture_output=True)
    common = [flow_cli, "--mesh", str(root/"mesh.solver.cm2d"), "--case", "channel", "--outlet-pressure", "101325",
              "--density", "1.176624281484062", "--pressure", "101325", "--u", "69.444", "--viscosity", "1.846e-5",
              "--conductivity", ".025759", "--wall-model", "no-slip", "--wall-thermal", "temperature", "--wall-value", "300",
              "--flux", "hllc", "--order", "2", "--integrator", "sdirk2", "--max-step", "1e-8", "--end-time", "1e-7", "--checkpoint-every", "1"]
    def run(name, options=(), expected=0):
        process = subprocess.run([*common, "--output", str(root/name), *options], capture_output=True, text=True, timeout=30)
        assert process.returncode == expected, process.stderr
        return json.loads((root/(name+".json")).read_text())
    tail = run("pressure-tail", ["--outlet-pressure", "101831.625", "--end-time", "7.200067830718531e-7"])
    assert tail["targetReached"] and tail["time"] == tail["requestedEndTime"]
    full = run("full")
    assert full["targetReached"] and not full["steadyConverged"]
    assert full["nonlinearTolerance"] == 2e-14
    # The optional stage budget must reach the native solver and survive an
    # accepted-state restart with the same control. This is not a PDE audit.
    algebraic = ["--nonlinear-tolerance", "1e-10"]
    candidate = run("algebraic", algebraic)
    assert candidate["targetReached"] and candidate["nonlinearTolerance"] == 1e-10
    assert 0 <= candidate["maximumAcceptedStageDefect"] <= 1e-10
    assert 0 <= candidate["maximumStageOutputDefect"] <= 8e-10
    run("algebraic-limited", [*algebraic, "--max-steps", "2"], 2)
    run("algebraic-resumed", [*algebraic, "--restart", str(root/"algebraic-limited.checkpoint")])
    assert (root/"algebraic.checkpoint").read_bytes() == (root/"algebraic-resumed.checkpoint").read_bytes()
    for index, options in enumerate((["--nonlinear-tolerance", value] for value in ["0", "-1", "1e-9", "nan", "inf"])):
        invalid_prefix = root/f"invalid-algebraic-{index}"
        rejected = subprocess.run([*common, "--output", str(invalid_prefix), *options], capture_output=True, text=True, timeout=30)
        assert rejected.returncode == 1 and not invalid_prefix.with_suffix(".json").exists()
    explicit_budget = subprocess.run([*common, "--output", str(root/"explicit-algebraic"), "--integrator", "explicit", *algebraic], capture_output=True, text=True, timeout=30)
    assert explicit_budget.returncode == 1 and "requires SDIRK2" in explicit_budget.stderr
    limited = run("limited", ["--max-steps", "2"], 2)
    assert limited["acceptedSteps"] == 2 and "budget" in limited["failure"]
    run("resumed", ["--restart", str(root/"limited.checkpoint")])
    assert (root/"full.checkpoint").read_bytes() == (root/"resumed.checkpoint").read_bytes()
    # Recovery of an unfinished desktop run must reuse the full native reader
    # without advancing or overwriting any result, even at/beyond the old goal.
    original_checkpoint = (root/"full.checkpoint").read_bytes()
    original_names = set(root.iterdir())
    def inspect(options=(), expected=0):
        process = subprocess.run([*common, "--inspect-restart", str(root/"full.checkpoint"), *options],
                                 cwd=root, capture_output=True, text=True, timeout=30)
        assert process.returncode == expected, process.stderr
        assert (root/"full.checkpoint").read_bytes() == original_checkpoint
        assert set(root.iterdir()) == original_names
        return process
    inspected = json.loads(inspect(["--end-time", "1e-10"]).stdout)
    assert inspected == {"type":"euler-checkpoint", "time":full["time"], "steps":full["steps"], "cells":full["cells"]}
    inspect(["--outlet-pressure", "101831.625"], 1)
    inspect(["--viscosity", "2e-5"], 1)
    inspect(["--output", str(root/"full")], 1)
    steady = run("steady-short", ["--mode", "steady", "--steady-scale", "5.76e-6", "--steady-tolerance", "1e-5"], 2)
    assert not steady["steadyConverged"] and not steady["targetReached"] and "horizon" in steady["failure"]
    budget = run("budget", ["--max-seconds", "0.000001"], 2)
    assert budget["steps"] == 0 and not budget["targetReached"]
    total = run("reservoir", ["--inlet-model", "total", "--inlet-total-pressure", "104190.2", "--inlet-total-temperature", "302.4"])
    assert total["targetReached"]
    run("reservoir-limited", ["--inlet-model", "total", "--inlet-total-pressure", "104190.2", "--inlet-total-temperature", "302.4", "--max-steps", "2"], 2)
    run("reservoir-resumed", ["--inlet-model", "total", "--inlet-total-pressure", "104190.2", "--inlet-total-temperature", "302.4", "--restart", str(root/"reservoir-limited.checkpoint")])
    assert (root/"reservoir.checkpoint").read_bytes() == (root/"reservoir-resumed.checkpoint").read_bytes()
    # Candidate/retry/restart semantics on a generated cut-cell channel with a
    # total inlet and a heated no-slip wall. Same-clock state agreement is
    # checked in the native euler_preconditioner test; this workflow may retry.
    heated = ["--inlet-model", "total", "--inlet-total-pressure", "104190.2", "--inlet-total-temperature", "302.4",
              "--wall-value", "330", "--max-step", "2e-8", "--end-time", "8e-8"]
    preconditioned = [*heated, "--implicit-preconditioner", "frozen-flux-ilu0"]
    ilu = run("heated-ilu", preconditioned)
    assert ilu["targetReached"] and ilu["implicitPreconditioner"] == "frozen-flux-ilu0"
    if ilu["rejectedCandidates"]:
        assert ilu["lastRejectedReason"]
    ilu_fields = json.loads((root/"heated-ilu.fields.json").read_text())["cells"]
    assert ilu["time"] == ilu["requestedEndTime"]
    assert len(ilu_fields) == ilu["cells"]
    assert all(c["rho"] > 0 and c["p"] > 0 and c["temperature"] > 0 for c in ilu_fields)
    run("heated-ilu-limited", [*preconditioned, "--max-steps", "2"], 2)
    run("heated-ilu-resumed", [*preconditioned, "--restart", str(root/"heated-ilu-limited.checkpoint")])
    assert (root/"heated-ilu.checkpoint").read_bytes() == (root/"heated-ilu-resumed.checkpoint").read_bytes()
    process = subprocess.Popen([*common, "--output", str(root/"cancelled"), "--end-time", "1e-4"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    line = process.stdout.readline()
    assert '"type":"euler-step"' in line
    process.send_signal(signal.SIGTERM)
    process.communicate(timeout=30)
    if sys.platform == "win32":
        # Windows TerminateProcess cannot call the C++ SIGTERM handler. Its
        # last atomic checkpoint survives; a final summary is not promised.
        assert process.returncode != 0
    else:
        assert process.returncode == 2
        cancelled = json.loads((root/"cancelled.json").read_text())
        assert cancelled["status"] == "cancelled" and cancelled["steps"] > 0
    inspection = subprocess.run([*common, "--inspect-restart", str(root/"cancelled.checkpoint")],
                                check=True, capture_output=True, text=True, timeout=30)
    saved = json.loads(inspection.stdout)
    assert saved["steps"] > 0 and saved["time"] > 0
    # The child may accept more steps before its parent is scheduled again.
    # Compare at one later physical endpoint, not an assumed cancellation step.
    end = saved["time"] + 1e-7
    run("cancel-full", ["--end-time", str(end)])
    run("cancel-resumed", ["--restart", str(root/"cancelled.checkpoint"), "--end-time", str(end)])
    assert (root/"cancel-full.checkpoint").read_bytes() == (root/"cancel-resumed.checkpoint").read_bytes()
print("native generated-grid cloud semantics passed")
