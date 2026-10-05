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
    limited = run("limited", ["--max-steps", "2"], 2)
    assert limited["acceptedSteps"] == 2 and "budget" in limited["failure"]
    run("resumed", ["--restart", str(root/"limited.checkpoint")])
    assert (root/"full.checkpoint").read_bytes() == (root/"resumed.checkpoint").read_bytes()
    steady = run("steady-short", ["--mode", "steady", "--steady-scale", "5.76e-6", "--steady-tolerance", "1e-5"], 2)
    assert not steady["steadyConverged"] and not steady["targetReached"] and "horizon" in steady["failure"]
    budget = run("budget", ["--max-seconds", "0.000001"], 2)
    assert budget["steps"] == 0 and not budget["targetReached"]
    total = run("reservoir", ["--inlet-model", "total", "--inlet-total-pressure", "104190.2", "--inlet-total-temperature", "302.4"])
    assert total["targetReached"]
    run("reservoir-limited", ["--inlet-model", "total", "--inlet-total-pressure", "104190.2", "--inlet-total-temperature", "302.4", "--max-steps", "2"], 2)
    run("reservoir-resumed", ["--inlet-model", "total", "--inlet-total-pressure", "104190.2", "--inlet-total-temperature", "302.4", "--restart", str(root/"reservoir-limited.checkpoint")])
    assert (root/"reservoir.checkpoint").read_bytes() == (root/"reservoir-resumed.checkpoint").read_bytes()
    process = subprocess.Popen([*common, "--output", str(root/"cancelled"), "--end-time", "1e-4"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    line = process.stdout.readline()
    assert '"type":"euler-step"' in line
    process.send_signal(signal.SIGTERM)
    process.communicate(timeout=30)
    assert process.returncode == 2
    cancelled = json.loads((root/"cancelled.json").read_text())
    assert cancelled["status"] == "cancelled" and cancelled["steps"] > 0
    run("cancel-resumed", ["--restart", str(root/"cancelled.checkpoint")])
    assert (root/"full.checkpoint").read_bytes() == (root/"cancel-resumed.checkpoint").read_bytes()
print("native generated-grid cloud semantics passed")
