"""External mature P2/P1 FEM reference on the frozen polygon.

Uses the analytic velocity trace family explicitly, not facet-wall motion.
This is an external solver experiment, not a duplicate native FV audit.
"""
from pathlib import Path
import argparse, importlib.metadata, json, resource, sys, time
p=argparse.ArgumentParser();p.add_argument('--mesh',type=Path,required=True)
p.add_argument('--output',type=Path,required=True);p.add_argument('--linear',choices=['direct','fieldsplit'],default='direct')
a=p.parse_args();sys.argv[:]=[sys.argv[0]];started=time.monotonic()
import firedrake as fd
import numpy as np

a.output.mkdir(parents=True,exist_ok=False)
record=dict(solver='Firedrake',version=importlib.metadata.version('firedrake'),
    boundaryModel='analytic-trace',element='Taylor-Hood CG2 velocity / CG1 pressure',
    mpiRanks=1,threads=1,gpuUsed=False,nu=.1,Uref=.5,linear=a.linear,
    scope='Independent mature FEM reference; differs from native cell-centred FV and facet-wall family')
def save():
    record.update(totalSeconds=time.monotonic()-started,
        cpuSeconds=resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime,
        peakRssMiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
    (a.output/'resource.json').write_text(json.dumps(record,indent=2)+'\n')
try:
    t=time.monotonic();mesh=fd.Mesh(str(a.mesh));V=fd.VectorFunctionSpace(mesh,'CG',2);Q=fd.FunctionSpace(mesh,'CG',1);W=V*Q
    x,y=fd.SpatialCoordinate(mesh);r2=x*x+y*y;factor=(1/r2-1)/3
    ue=fd.as_vector((-y*factor,x*factor));pe=(r2-1/r2)/18-fd.ln(r2)/9
    solution=fd.Function(W);u,pressure=fd.split(solution);v,q=fd.TestFunctions(W)
    form=(.1*fd.inner(fd.grad(u)+fd.grad(u).T,fd.grad(v))+fd.inner(fd.dot(fd.grad(u),u),v)
          -pressure*fd.div(v)+q*fd.div(u))*fd.dx
    boundary=fd.Function(V).interpolate(ue)
    bcs=[fd.DirichletBC(W.sub(0),boundary,'on_boundary')]
    null=fd.MixedVectorSpaceBasis(W,[W.sub(0),fd.VectorSpaceBasis(constant=True)])
    options={'snes_type':'newtonls','snes_linesearch_type':'bt','snes_rtol':1e-8,'snes_atol':1e-10,
        'snes_max_it':25,'snes_monitor':None,'ksp_monitor':None,
        'ksp_type':'preonly','pc_type':'lu','pc_factor_mat_solver_type':'mumps',
        'mat_mumps_icntl_14':80}
    if a.linear=='fieldsplit':
        options.update(ksp_type='fgmres',ksp_rtol=1e-9,pc_type='fieldsplit',pc_fieldsplit_type='schur',
            pc_fieldsplit_schur_fact_type='full',fieldsplit_0_ksp_type='preonly',fieldsplit_0_pc_type='hypre',
            fieldsplit_1_ksp_type='preonly',fieldsplit_1_pc_type='jacobi')
        options.pop('pc_factor_mat_solver_type')
    solver=fd.NonlinearVariationalSolver(fd.NonlinearVariationalProblem(form,solution,bcs=bcs),
        solver_parameters=options,nullspace=null,transpose_nullspace=null)
    record.update(cells=mesh.num_cells(),velocityDofs=V.dim(),pressureDofs=Q.dim(),totalDofs=W.dim(),
        setupSeconds=time.monotonic()-t,solverParameters=options);save()
    t=time.monotonic();solver.solve();record['solveSeconds']=time.monotonic()-t
    uh,ph=solution.subfunctions;area=fd.assemble(fd.Constant(1)*fd.dx(domain=mesh))
    gauge=fd.assemble((ph-pe)*fd.dx)/area
    span=.0543006265422344
    record.update(converged=solver.snes.getConvergedReason()>0,snesReason=solver.snes.getConvergedReason(),
        nonlinearIterations=solver.snes.getIterationNumber(),linearIterations=solver.snes.getLinearSolveIterations(),
        nonlinearResidual=solver.snes.getFunctionNorm(),
        velocityRmsOverUref=float(fd.sqrt(fd.assemble(fd.inner(uh-ue,uh-ue)*fd.dx)/area)/.5),
        pressureRmsOverAnalyticSpan=float(fd.sqrt(fd.assemble((ph-pe-gauge)**2*fd.dx)/area)/span),
        netBoundaryFlux=float(fd.assemble(fd.dot(uh,fd.FacetNormal(mesh))*fd.ds)),area=float(area))
    # Keep one restartable mature-solver state rather than thousands of snapshots.
    with fd.CheckpointFile(str(a.output/'accepted.h5'),'w') as f:
        f.save_mesh(mesh);f.save_function(uh,name='velocity');f.save_function(ph,name='pressure')
    record['status']='finished'
except Exception as error:
    record.update(status='failed',converged=False,error=repr(error));save();raise
finally:save()
print(json.dumps(record))
