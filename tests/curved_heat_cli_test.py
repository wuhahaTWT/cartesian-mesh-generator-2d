"""Native real-grid annular spatial regression, not an independent PDE audit."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

repo, build = map(Path, sys.argv[1:3])
with tempfile.TemporaryDirectory(prefix="cartmesh-curved-") as d:
    output = Path(d)/"run"
    subprocess.run([sys.executable, str(repo/"tools/flow/run_curved_wall.py"), "--build", str(build),
                    "--output", str(output), "--skip-coupled"], check=True, timeout=60)
    j = json.loads((output/"run.json").read_text())
    a, b, c = [j["cases"][f"L{level}-linear-trace"] for level in [4, 5, 6]]
    # This specific 128-facet annulus and requested resolutions: <2% wall heat
    # L1, normalized by exact absolute chord heat, and improvement over coarse.
    # Not a universal curved-wall physical accuracy standard.
    assert c["wallFluxRelativeL1"] < .02 and c["wallFluxRelativeL1"] < a["wallFluxRelativeL1"]/2
    assert c["temperatureRelativeL1"] < b["temperatureRelativeL1"] < a["temperatureRelativeL1"]
    for x in [a, b, c]:
        assert x["nonpositiveDiffusionDiagonalRows"] == 0
        assert x["algebraicResidualRelative"] < 1e-7*x["wallFluxRelativeL1"]
        assert x["toleranceSensitivityK"] < .01*.2*x["temperatureRelativeL1"]
        assert x["innerHeat"] > 0 and x["outerHeat"] < 0
print("native annular linear-wall spatial regression passed; quadratic qualification excluded")
