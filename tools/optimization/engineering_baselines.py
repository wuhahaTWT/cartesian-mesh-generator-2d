"""Parametric engineering reference families, before any topology optimization.

These generate centreline geometry, not a smoothed version of an optimized XY.
Widths are fitted to the same filtered/projected fluid-area constraint. Every
family member is evaluated; the lowest objective is the primary reference.
"""
import numpy as np


def distance_to_path(points, path):
    distance = np.full(len(points), np.inf)
    for a,b in zip(path[:-1],path[1:]):
        ab=b-a
        t=np.clip((points-a)@ab/(ab@ab),0,1)
        distance=np.minimum(distance,np.linalg.norm(points-a-t[:,None]*ab,axis=1))
    return distance


def paths(model, parameter):
    p=model.problem; w,h=p.width,p.height
    if p.case=="elbow":
        r=parameter*min(w,h)
        c=np.array([.8*w-r,.8*h-r])
        angle=np.linspace(np.pi/2,0,65)
        curve=c+r*np.c_[np.cos(angle),np.sin(angle)]
        return [np.vstack(([0,.8*h],curve,[.8*w,0]))]
    if p.case=="four-terminal":
        # Explicit U-return or parallel-channel alternatives. Neither is called
        # a topology-optimized design, and both are retained in the report.
        if parameter==0:
            return [np.array([[0,y*h],[w,y*h]]) for y in (.3,.7)]
        reach=parameter*w
        angle=np.linspace(np.pi/2,-np.pi/2,65)
        arc=np.c_[reach*np.cos(angle),.5*h+.2*h*np.sin(angle)]
        return [arc, np.c_[w-arc[:,0],arc[:,1]]]
    if p.case=="double-pipe":
        t=parameter*w
        return [np.array([[0,y*h],[t,.5*h],[w-t,.5*h],[w,y*h]]) for y in (.25,.75)]
    if p.case=="diffuser":
        return [np.array([[0,.5*h],[w,.5*h]])]
    return [np.array([[0,.25*h],[w,.75*h]])]


def family(model,beta):
    p=model.problem
    xx,yy=np.meshgrid((np.arange(model.nx)+.5)*model.dx,(np.arange(model.ny)+.5)*model.dy)
    points=np.c_[xx.ravel(),yy.ravel()]
    parameters={"elbow":[.3,.5,.7],"double-pipe":[.25,.35,.45],
                "diffuser":[.65,1.,1.5],"four-terminal":[0,.25,.4]}.get(p.case,[1.])
    result=[]
    for parameter in parameters:
        if p.case=="diffuser":
            # Symmetric expanding duct, three power-law taper rates.
            t=points[:,0]/p.width
            base_half=(p.port_width+(p.height-p.port_width)*t**parameter)/2
            signed=abs(points[:,1]-p.height/2)-base_half
        else:
            signed=np.minimum.reduce([distance_to_path(points,path) for path in paths(model,parameter)])
        def candidate(offset):
            return model.enforce_passive(np.clip(.5+(offset-signed)/min(model.dx,model.dy),0,1))
        lo,hi=-max(p.width,p.height),max(p.width,p.height)
        for _ in range(55):
            mid=(lo+hi)/2
            if model.volume(candidate(mid),beta)[0]>p.volume_fraction:
                hi=mid
            else:
                lo=mid
        x=candidate(lo)
        if abs(model.volume(x,beta)[0]-p.volume_fraction)>1e-10:
            raise ValueError("engineering baseline volume infeasible")
        result.append((f"{p.case}-parameter-{parameter:g}",x))
    return result


def best_reference(model,beta,q=.1,objective="pressure-power"):
    rows=[]; best=None
    for label,x in family(model,beta):
        try:
            e=model.evaluate(x,q,beta,objective)
            rows.append(dict(label=label,objective=e.objective,volume=e.volume,status="analyzed"))
            if best is None or e.objective<best[2].objective:
                best=(label,x,e)
        except (ArithmeticError,ValueError,RuntimeError) as exc:
            rows.append(dict(label=label,status="analysis-failed",issue=str(exc)))
    if best is None:
        raise ArithmeticError("all engineering reference analyses failed")
    return best,rows
