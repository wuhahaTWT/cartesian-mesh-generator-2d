"""Transparent leakage, three-grid GCI and ranking statistics for T01."""
import math
import numpy as np
from scipy.optimize import brentq
from scipy.stats import kendalltau


def leakage(model,e,q):
    rho=e.rho
    # Dual-face absolute volume throughput, counted once per internal face.
    # It is a sampling measure, NOT the fraction of distinct inlet particles
    # that penetrate a solid; recirculation can cross several sampling faces.
    solid=model.drag_weights@rho/np.asarray(model.drag_weights.sum(axis=1)).ravel()<.5
    lengths=np.r_[np.full(model.nu_faces,model.dy),np.full(model.nv_faces,model.dx)]
    flux=np.abs(e.velocity[model.free])*lengths[model.free]
    fraction=float(np.sum(flux[solid[model.free]])/max(np.sum(flux),1e-30))
    alpha=model.problem.alpha_max*q*(1-rho)/(q+rho)
    drag=alpha*np.asarray(model.drag_weights.T@(e.velocity**2)).ravel()
    k=model.viscous.tocoo()
    off=k.row!=k.col
    face_energy=np.bincount(k.row[off],weights=-.5*k.data[off]*(e.velocity[k.row[off]]-e.velocity[k.col[off]])**2,
                            minlength=model.velocities)
    face_energy+=np.asarray(model.viscous.sum(axis=1)).ravel()*e.velocity**2
    dual_area=np.asarray(model.drag_weights.sum(axis=1)).ravel()
    viscous=np.asarray(model.drag_weights.T@(face_energy/dual_area)).ravel()
    dissipation=drag+viscous
    return dict(solidFluxFraction=fraction,
                solidDissipationFraction=float(np.sum(dissipation[rho<.5])/max(e.dissipation,1e-30)),
                solidDragDissipation=float(np.sum(drag[rho<.5])),
                dissipationReconstructionRelativeError=float(abs(np.sum(dissipation)-e.dissipation)/max(abs(e.dissipation),1e-30)),
                definition="Absolute internal MAC face volume throughput through mean face rho<0.5 / all absolute internal face throughput; solid-cell share of positive drag+viscous energy; repeated crossings counted.")


def grid_uncertainty(rows):
    accepted=[r for r in rows if r.get("metrics") is not None]
    accepted=sorted(accepted,key=lambda r:r["level"])
    result=dict(method="unavailable",acceptedLevels=[r["level"] for r in accepted],uncertainty=None,
                allAttemptedGridsAccepted=len(accepted)==len(rows),physicalAccuracyQualified=False)
    if not accepted:return result
    result["J_T"]=accepted[-1]["metrics"]["dimensionlessTotalPower"]
    if len(accepted)<2:return result
    values=[r["metrics"]["dimensionlessTotalPower"] for r in accepted]
    difference=abs(values[-1]-values[-2])
    result.update(method="two-grid-change",twoGridChange=difference,
                  twoGridRelativeChange=difference/max(abs(values[-1]),1e-30))
    if len(accepted)<3:return result
    coarse,medium,fine=accepted[-3:]
    families=[(r.get('meshPaddingFraction',1/30),r.get('meshAggregationFraction',.25)) for r in (coarse,medium,fine)]
    if len(set(families))!=1:
        result['gciUnavailableReason']='mixed background-grid padding families or aggregation settings; retain observed grid changes without a three-grid GCI fit'
        return result
    h3,h2,h1=[math.sqrt(r["metrics"]["area"]/r["metrics"]["cells"]) for r in (coarse,medium,fine)]
    r21,r32=h2/h1,h3/h2
    d32,d21=values[-3]-values[-2],values[-2]-values[-1]
    result["effectiveGridSpacings"]=[h3,h2,h1]
    result["refinementRatios"]=[r32,r21]
    if r21<=1 or r32<=1 or d32*d21<=0:
        result["gciUnavailableReason"]="nonmonotone/zero differences or nonrefining effective spacing"
        return result
    ratio=abs(d32/d21)
    def equation(p):
        return r21**p*math.expm1(p*math.log(r32))/math.expm1(p*math.log(r21))-ratio
    try:
        order=brentq(equation,1e-6,64)
    except (ValueError,OverflowError):
        result["gciUnavailableReason"]="no positive observed-order root for the three grids"
        return result
    uncertainty=1.25*abs(d21)/math.expm1(order*math.log(r21))
    result.update(method="three-grid-GCI",observedOrder=order,safetyFactor=1.25,
                  uncertainty=uncertainty,gciRelative=uncertainty/max(abs(values[-1]),1e-30),
                  extrapolated=values[-1]+(values[-1]-values[-2])/math.expm1(order*math.log(r21)),
                  asymptoticRangeIndependentlyEstablished=False,
                  caveat="Roache/Richardson estimate conditional on an asymptotic power law; three data alone do not establish that law. No physical-accuracy qualification.")
    return result


def paired_grid_orderings(a,b):
    """Compare only shared native levels; stability is not an error bound."""
    def values(row):
        return {(r['level'],r.get('meshPaddingFraction',1/30),r.get('meshAggregationFraction',.25)):r['metrics']['dimensionlessTotalPower'] for r in row.get('native',[])
                if r.get('metrics') is not None}
    av,bv=values(a),values(b)
    db=a['porousObjective']-b['porousObjective']
    shared=sorted(av.keys()&bv.keys())
    comparisons=[dict(level=l[0],meshPaddingFraction=l[1],meshAggregationFraction=l[2],nativeA=av[l],nativeB=bv[l],difference=av[l]-bv[l],
                      reversed=db*(av[l]-bv[l])<0) for l in shared]
    return dict(commonLevels=[l[0] for l in shared],commonGridContracts=[dict(level=l[0],meshPaddingFraction=l[1],meshAggregationFraction=l[2]) for l in shared],comparisons=comparisons,
                reversalObservedOnAtLeastTwoCommonGrids=len(comparisons)>=2 and all(x['reversed'] for x in comparisons),
                nativeOrderChangesAcrossCommonGrids=len({np.sign(x['difference']) for x in comparisons})>1,
                interpretation='Observed ordering on identical native levels; neither grid independence nor a confidence interval.')


def ranking(rows):
    available=[r for r in rows if r.get("grid",{}).get("J_T") is not None and r.get("porousObjective") is not None]
    qualified=[r for r in available if r["grid"].get("uncertainty") is not None]
    result=dict(attempted=len(rows),available=len(available),gciAvailable=len(qualified),
                kendallTau=None,allPairs=[],reversals=[],conditionalGciReversals=[],
                resolvedReversals=[],twoGridSensitivityReversals=[],
                interpretation="Conditional GCI comparisons are observations. Resolved rankings require independently established native grid uncertainty and one porous evaluation contract.")
    if len(available)>=2:
        tau=kendalltau([r["porousObjective"] for r in available],[r["grid"]["J_T"] for r in available]).statistic
        result["kendallTau"]=float(tau) if np.isfinite(tau) else None
    for name,subset in (("atLeastTwoGrids",[r for r in available if len(r["grid"].get("acceptedLevels",[]))>=2]),
                        ("threeGridGci",qualified)):
        tau=kendalltau([r["porousObjective"] for r in subset],[r["grid"]["J_T"] for r in subset]).statistic if len(subset)>=2 else None
        result[name]=dict(count=len(subset),kendallTau=float(tau) if tau is not None and np.isfinite(tau) else None)
    for i,a in enumerate(available):
        for b in available[i+1:]:
            db=a["porousObjective"]-b["porousObjective"]
            dt=a["grid"]["J_T"]-b["grid"]["J_T"]
            pair=dict(a=a["id"],b=b["id"],porousDifference=db,trueDifference=dt,reversed=db*dt<0)
            ua,ub=a["grid"].get("uncertainty"),b["grid"].get("uncertainty")
            pair["uncertaintySum"]=ua+ub if ua is not None and ub is not None else None
            pair["exceedsConditionalGci"]=pair["reversed"] and pair["uncertaintySum"] is not None and abs(dt)>pair["uncertaintySum"]
            def contract(row):
                if "porousEvaluationContract" in row:
                    return row["porousEvaluationContract"]
                problem=row.get("problem",{})
                keys=("case","width","height","port_width","viscosity","reynolds","alpha_max","nx","ny")
                final=row.get("finalParameters",{})
                if not all(k in problem for k in keys) or not all(k in final for k in ("q","beta")):
                    return None
                return dict(problem={k:problem[k] for k in keys},q=final["q"],beta=final["beta"])
            ac,bc=contract(a),contract(b)
            pair["commonPorousModelEstablished"]=ac is not None and bc is not None and ac==bc
            pair["nativeUncertaintyEstablished"]=all(r["grid"].get("asymptoticRangeIndependentlyEstablished") is True for r in (a,b))
            pair["resolved"]=bool(pair["exceedsConditionalGci"] and pair["commonPorousModelEstablished"] and pair["nativeUncertaintyEstablished"])
            pair['pairedNativeGrids']=paired_grid_orderings(a,b)
            result["allPairs"].append(pair)
            if pair["reversed"]:result["reversals"].append(pair)
            if pair["exceedsConditionalGci"]:result["conditionalGciReversals"].append(pair)
            if pair["resolved"]:result["resolvedReversals"].append(pair)
            ca,cb=a["grid"].get("twoGridChange"),b["grid"].get("twoGridChange")
            if pair["reversed"] and ca is not None and cb is not None and abs(dt)>ca+cb:
                result["twoGridSensitivityReversals"].append(pair)
    return result
