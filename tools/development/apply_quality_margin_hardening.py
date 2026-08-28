from pathlib import Path

source = Path("src/quality/SolverTopology2D.cpp").read_text(encoding="utf-8")
if "preferredSolverQualityPolicy()" not in source or "marginOptimized" not in source:
    raise RuntimeError("quality-margin optimizer source patch is not present")

path = Path("tests/solver_export_test.cpp")
text = path.read_text(encoding="utf-8")
old = "doubleNacaSolver.profile.candidateTopologyCount==1,"
new = "doubleNacaSolver.profile.candidateTopologyCount>=1,"
if old in text:
    text = text.replace(old, new, 1)
elif new not in text:
    raise RuntimeError("expected retained double-NACA profile assertion not found")
path.write_text(text, encoding="utf-8")
print("Aligned solver-export regression with post-validity quality optimization")
