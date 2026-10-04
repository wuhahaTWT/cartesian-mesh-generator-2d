"""Explicit benchmark port contracts and literature mappings for T01.

All boundary profiles are prescribed on both inlet and outlet in the MAC model.
Literature-inspired variants are labelled as adaptations, never exact replicas.
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class Port:
    side: str
    centre: float
    width: float
    peak: float
    role: str
    name: str


def ports_for(p):
    get = (lambda key: p[key]) if isinstance(p, dict) else (lambda key: getattr(p, key))
    case, h, w = get("case"), get("height"), get("port_width")
    if case == "channel":
        return [Port("left",h/2,h,1,"inlet","inlet_0"), Port("right",h/2,h,1,"outlet","outlet_0")]
    if case == "diffuser":
        return [Port("left",h/2,w,1,"inlet","inlet_0"), Port("right",h/2,h,w/h,"outlet","outlet_0")]
    if case == "elbow":
        return [Port("left",.8*h,w,1,"inlet","inlet_0"), Port("bottom",.8*get("width"),w,1,"outlet","outlet_0")]
    if case == "four-terminal":
        return [Port("left",.7*h,w,1,"inlet","inlet_0"), Port("right",.3*h,w,1,"inlet","inlet_1"),
                Port("left",.3*h,w,1,"outlet","outlet_0"), Port("right",.7*h,w,1,"outlet","outlet_1")]
    left = [h/4,3*h/4] if case == "double-pipe" else [h/4]
    right = left if case == "double-pipe" else [3*h/4]
    return [Port("left",c,w,1,"inlet",f"inlet_{i}") for i,c in enumerate(left)]+[
            Port("right",c,w,1,"outlet",f"outlet_{i}") for i,c in enumerate(right)]


CASES = {
    "double-pipe": dict(width=1.5, height=1., port_width=1/6, volume_fraction=1/3, nx=36, ny=24,
        source="https://doi.org/10.1002/fld.426",
        accessibleParameters="https://www.dolfin-adjoint.org/en/stable/documentation/stokes-topology/stokes-topology.html",
        mapping="Borrvall-Petersson Fig.10: H=1, L=1.5, w=1/6, centres 1/4 and 3/4, volume 1/3. MAC/filter/passive collars differ from FEM.",
        referenceTopology="A common trunk can outperform two separate paths at this aspect ratio.",
        quantitativeReference=None),
    "diffuser": dict(width=1., height=1., port_width=1/3, volume_fraction=.5, nx=24, ny=24,
        source="https://doi.org/10.1002/fld.426",
        mapping="Diffuser-inspired square with inlet w=1/3, outlet H=1 and equal flow (outlet peak=1/3), volume 1/2. Original full-text parameters unavailable; this is explicitly an adapted benchmark, not a numerical reproduction.",
        referenceTopology="Smooth expansion is the engineering reference; exact published contour/numbers not verified.",
        quantitativeReference=None),
    "elbow": dict(width=1., height=1., port_width=.2, volume_fraction=.25, nx=30, ny=30,
        source="https://doi.org/10.1002/fld.426",
        inertiaSource="https://doi.org/10.1007/s00158-004-0508-7",
        accessibleParameters="https://backend.orbit.dtu.dk/ws/portalfiles/portal/4897055/C_Agh_PhD_Texts_Thesis_070122b_final_070410_final_Forside_AGH_thesis.pdf",
        mapping="Borrvall-Petersson bend as drawn in Gersborg-Hansen thesis Fig.C.1: divide coordinates by 5; left inlet y=0.8, bottom outlet x=0.8, width 0.2. Volume 0.25. Zero fluid drag differs from the finite fluid drag in the 2005 inertia study.",
        referenceTopology="Low-Re diagonal shortcut; curvature may increase with inertia. No exact objective reference under these controls.",
        quantitativeReference=None),
    "four-terminal": dict(width=.7, height=1., port_width=.2, volume_fraction=.4, nx=21, ny=30,
        source="https://doi.org/10.1002/nme.1468", accessibleParameters="https://arxiv.org/pdf/physics/0410086",
        mapping="Olesen et al. Fig.6-7: H=5 ell, L=3.5 ell -> H=1,L=0.7,w=0.2, centres 0.3/0.7, volume 0.4. Opposite diagonal inlets; two U-turns compete with straight paths. Adaptation omits 2 ell leads and prescribes outlet profiles instead of pressure. Re_mean=(2/3)Re_peak, so paper 20/200 map to 13.333/133.333.",
        referenceTopology="Published low-Re two U-turns versus high-Re parallel paths; exact crossover is not transferable after the boundary adaptation.",
        quantitativeReference=None),
}


def problem_parameters(case):
    return {k:v for k,v in CASES[case].items() if k in
            ("width","height","port_width","volume_fraction","nx","ny")}
