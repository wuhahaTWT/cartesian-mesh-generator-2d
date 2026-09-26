"""Experimental native-2D staggered Stokes--Brinkman topology analysis.

This is a dedicated optimisation discretisation, NOT the product background-grid
format and NOT a replacement for the polygonal incompressible solver. Variables
are nondimensional. rho=1 is fluid; rho=0 is a finite-resistance approximation to
solid. Extracted walls require a separate sharp-interface CFD calculation.

The symmetric, Dirichlet-eliminated system is [K, -D.T; -D, 0]. K contains
viscosity and diagonal Brinkman resistance; D is the integrated MAC divergence.
One redundant continuity equation / pressure constant is removed. Gradients use
the transpose of THIS discrete system, including the filter and projection.
"""
from dataclasses import dataclass
import math

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu


@dataclass(frozen=True)
class Problem:
    nx: int = 48
    ny: int = 32
    width: float = 1.5
    height: float = 1.0
    viscosity: float = 1.0
    alpha_max: float = 25000.0
    volume_fraction: float = 1.0 / 3.0
    filter_radius: float = 0.06
    port_width: float = 1.0 / 6.0
    case: str = "double-pipe"
    # Algebraic residual / rhs infinity norm, not a physical accuracy claim.
    linear_residual_limit: float = 1e-8

    def __post_init__(self):
        if (not isinstance(self.nx, int) or not isinstance(self.ny, int) or
                not 8 <= self.nx <= 160 or not 8 <= self.ny <= 160):
            raise ValueError("nx and ny must be integers in [8,160]")
        for name in ("width", "height", "viscosity", "filter_radius", "port_width",
                     "linear_residual_limit"):
            if not math.isfinite(getattr(self, name)) or getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive and finite")
        if not math.isfinite(self.alpha_max) or self.alpha_max < 0:
            raise ValueError("alpha_max must be nonnegative and finite")
        if not 0 < self.volume_fraction <= 1:
            raise ValueError("volume_fraction must be in (0,1]")
        if self.case not in ("double-pipe", "bend", "channel"):
            raise ValueError("unknown flow topology case")
        if self.case != "channel" and not 2*self.height/self.ny <= self.port_width < self.height/2:
            raise ValueError("ports need at least two cells and must be narrower than half-height")
        if self.filter_radius > min(self.width, self.height)/4:
            raise ValueError("filter radius exceeds the bounded prototype range")


@dataclass
class Evaluation:
    objective: float
    gradient: np.ndarray
    volume: float
    volume_gradient: np.ndarray
    rho: np.ndarray
    velocity: np.ndarray
    pressure: np.ndarray
    linear_residual: float
    continuity: float
    adjoint_residual: float
    dissipation: float
    pressure_power: float
    grayness: float
    solid_speed_fraction: float


def port_average(y0, y1, centre, width):
    """Exact face average of a unit-peak parabolic aperture velocity."""
    lo, hi = max(y0, centre-width/2), min(y1, centre+width/2)
    if hi <= lo:
        return 0.0
    primitive = lambda y: y - 4*(y-centre)**3/(3*width**2)
    return (primitive(hi)-primitive(lo))/(y1-y0)


class StokesBrinkman:
    def __init__(self, problem=Problem()):
        self.problem = problem
        p = problem
        self.nx, self.ny = p.nx, p.ny
        self.dx, self.dy = p.width/p.nx, p.height/p.ny
        self.area = self.dx*self.dy
        self.cells = p.nx*p.ny
        self.nu_faces = (p.nx+1)*p.ny
        self.nv_faces = p.nx*(p.ny+1)
        self.velocities = self.nu_faces+self.nv_faces
        self.fixed_velocity = np.zeros(self.velocities)
        self.is_dirichlet = np.zeros(self.velocities, dtype=bool)
        self.fixed_design = np.full(self.cells, np.nan)
        self._boundary_data()
        self.free = np.flatnonzero(~self.is_dirichlet)
        self.prescribed = np.flatnonzero(self.is_dirichlet)
        self.design = np.isnan(self.fixed_design)
        self._assemble()
        self.filter = self._density_filter()
        self.pressure_weights = np.zeros(self.cells)
        for j in range(p.ny):
            self.pressure_weights[self.c(0, j)] += self.fixed_velocity[self.u(0, j)]*self.dy
            self.pressure_weights[self.c(p.nx-1, j)] -= self.fixed_velocity[self.u(p.nx, j)]*self.dy
        self.inflow = sum(self.fixed_velocity[self.u(0, j)]*self.dy for j in range(p.ny))
        outflow = sum(self.fixed_velocity[self.u(p.nx, j)]*self.dy for j in range(p.ny))
        if self.inflow <= 0 or abs(self.inflow-outflow) > 1e-12*self.inflow:
            raise ValueError("prescribed inlet and outlet fluxes must balance")

    def c(self, i, j):
        return j*self.nx+i

    def u(self, i, j):
        return j*(self.nx+1)+i

    def v(self, i, j):
        return self.nu_faces+j*self.nx+i

    def _boundary_data(self):
        p = self.problem
        left = [p.height/4, 3*p.height/4] if p.case == "double-pipe" else [p.height/4]
        right = left if p.case == "double-pipe" else [3*p.height/4]
        width = p.port_width
        if p.case == "channel":
            left = right = [p.height/2]
            width = p.height
        self.left_ports, self.right_ports, self.aperture_width = left, right, width
        for j in range(self.ny):
            for i, ports in ((0, left), (self.nx, right)):
                k = self.u(i, j)
                self.is_dirichlet[k] = True
                self.fixed_velocity[k] = sum(port_average(j*self.dy, (j+1)*self.dy, y, width)
                                             for y in ports)
        for i in range(self.nx):
            self.is_dirichlet[self.v(i, 0)] = self.is_dirichlet[self.v(i, self.ny)] = True
        if p.case == "channel":
            return
        # A two-cell passive collar preserves the physical openings. Wall cells
        # are fixed solid. These are geometric boundary conditions, not a seed
        # of the interior flow path, whose connectivity remains unconstrained.
        for j in range(self.ny):
            for i in range(self.nx):
                k = self.c(i, j)
                if j in (0, self.ny-1):
                    self.fixed_design[k] = 0.0
                if i < 2 or i >= self.nx-2:
                    side = 0 if i < 2 else self.nx
                    self.fixed_design[k] = float(self.fixed_velocity[self.u(side, j)] > 0)

    def _assemble(self):
        row, col, val = [], [], []
        def add(i, j, a):
            row.append(i); col.append(j); val.append(a)
        def spring(a, b, conductance):
            add(a, a, conductance); add(b, b, conductance)
            add(a, b, -conductance); add(b, a, -conductance)
        mu, dx, dy = self.problem.viscosity, self.dx, self.dy
        for j in range(self.ny):
            for i in range(self.nx):
                spring(self.u(i, j), self.u(i+1, j), mu*dy/dx)
        for i in range(self.nx+1):
            dual_width = dx*(0.5 if i in (0, self.nx) else 1)
            for j in range(self.ny-1):
                spring(self.u(i, j), self.u(i, j+1), mu*dual_width/dy)
            add(self.u(i, 0), self.u(i, 0), 2*mu*dual_width/dy)
            add(self.u(i, self.ny-1), self.u(i, self.ny-1), 2*mu*dual_width/dy)
        for i in range(self.nx):
            for j in range(self.ny):
                spring(self.v(i, j), self.v(i, j+1), mu*dx/dy)
        for j in range(self.ny+1):
            dual_height = dy*(0.5 if j in (0, self.ny) else 1)
            for i in range(self.nx-1):
                spring(self.v(i, j), self.v(i+1, j), mu*dual_height/dx)
            add(self.v(0, j), self.v(0, j), 2*mu*dual_height/dx)
            add(self.v(self.nx-1, j), self.v(self.nx-1, j), 2*mu*dual_height/dx)
        self.viscous = sparse.coo_matrix((val, (row, col)), shape=(self.velocities,)*2).tocsc()
        dr, dc, dv, wr, wc, wv = [], [], [], [], [], []
        for j in range(self.ny):
            for i in range(self.nx):
                k = self.c(i, j)
                for face, sign, length in ((self.u(i, j), -1, dy), (self.u(i+1, j), 1, dy),
                                            (self.v(i, j), -1, dx), (self.v(i, j+1), 1, dx)):
                    dr.append(k); dc.append(face); dv.append(sign*length)
                    wr.append(face); wc.append(k); wv.append(self.area/2)
        self.divergence = sparse.coo_matrix((dv, (dr, dc)), shape=(self.cells, self.velocities)).tocsr()
        self.drag_weights = sparse.coo_matrix((wv, (wr, wc)), shape=(self.velocities, self.cells)).tocsr()
        self.free_weights = self.drag_weights[self.free, :]
        self.dfree = self.divergence[:-1, self.free].tocsc()
        self.rhs = np.r_[-(self.viscous[self.free, :]@self.fixed_velocity),
                         self.divergence[:-1, :]@self.fixed_velocity]

    def _density_filter(self):
        radius = self.problem.filter_radius
        ir, jr = math.ceil(radius/self.dx), math.ceil(radius/self.dy)
        row, col, val = [], [], []
        offsets = [(di, dj, radius-math.hypot(di*self.dx, dj*self.dy))
                   for dj in range(-jr, jr+1) for di in range(-ir, ir+1)
                   if math.hypot(di*self.dx, dj*self.dy) < radius]
        for j in range(self.ny):
            for i in range(self.nx):
                for di, dj, w in offsets:
                    if 0 <= i+di < self.nx and 0 <= j+dj < self.ny:
                        row.append(self.c(i, j)); col.append(self.c(i+di, j+dj)); val.append(w)
        h = sparse.coo_matrix((val, (row, col)), shape=(self.cells,)*2).tocsr()
        return sparse.diags(1/np.asarray(h.sum(axis=1)).ravel())@h

    def enforce_passive(self, x):
        x = np.array(x, dtype=float, copy=True)
        if x.shape != (self.cells,) or not np.isfinite(x).all():
            raise ValueError("design must contain one finite value per analysis cell")
        x[~self.design] = self.fixed_design[~self.design]
        if np.any(x < 0) or np.any(x > 1):
            raise ValueError("design outside [0,1]")
        return x

    def physical(self, x, beta):
        x = self.enforce_passive(x)
        if not math.isfinite(beta) or not 0 <= beta <= 32:
            raise ValueError("projection beta must be in [0,32]")
        filtered = self.filter@x
        if beta == 0:
            rho, derivative = filtered.copy(), np.ones(self.cells)
        else:
            denominator = 2*np.tanh(beta/2)
            t = np.tanh(beta*(filtered-0.5))
            rho = (np.tanh(beta/2)+t)/denominator
            derivative = beta*(1-t*t)/denominator
        rho[~self.design] = self.fixed_design[~self.design]
        derivative[~self.design] = 0
        return rho, derivative

    def chain(self, derivative, physical_gradient):
        result = np.asarray(self.filter.T@(derivative*physical_gradient)).ravel()
        result[~self.design] = 0
        return result

    def volume(self, x, beta):
        rho, derivative = self.physical(x, beta)
        return float(np.mean(rho)), self.chain(derivative, np.full(self.cells, 1/self.cells))

    def feasible_design(self, x, beta):
        """Restore the physical (filtered/projected) volume by scalar bisection."""
        x = self.enforce_passive(x)
        target = self.problem.volume_fraction
        volume, _ = self.volume(x, beta)
        if abs(volume-target) <= 1e-12:
            return x
        end = x.copy()
        end[self.design] = float(volume < target)
        endpoint, _ = self.volume(end, beta)
        if min(volume, endpoint) > target+1e-12 or max(volume, endpoint) < target-1e-12:
            raise ValueError("volume fraction is infeasible with passive collars and filter")
        lo, hi = 0.0, 1.0
        for _ in range(50):
            mid = (lo+hi)/2
            candidate = (1-mid)*x+mid*end
            value, _ = self.volume(candidate, beta)
            if (value < target) == (volume < target):
                lo = mid
            else:
                hi = mid
        # Choose the feasible side, including for continuation-stage changes.
        weight = lo if volume < target else hi
        return self.enforce_passive((1-weight)*x+weight*end)

    def uniform_design(self, beta=0):
        return self.feasible_design(np.full(self.cells, self.problem.volume_fraction), beta)

    def evaluate(self, x, q=0.1, beta=0.0, objective="dissipation"):
        if not math.isfinite(q) or q <= 0:
            raise ValueError("Brinkman interpolation q must be positive")
        if objective not in ("dissipation", "pressure-power"):
            raise ValueError("unknown objective")
        rho, projection_derivative = self.physical(x, beta)
        alpha = self.problem.alpha_max*q*(1-rho)/(q+rho)
        dalpha = -self.problem.alpha_max*q*(1+q)/(q+rho)**2
        resistance = self.drag_weights@alpha
        stiffness = self.viscous+sparse.diags(resistance, format="csc")
        kii = stiffness[self.free, :][:, self.free]
        matrix = sparse.bmat([[kii, -self.dfree.T], [-self.dfree, None]], format="csc")
        factor = splu(matrix)
        state = factor.solve(self.rhs)
        linear = float(np.max(np.abs(matrix@state-self.rhs))/max(np.max(np.abs(self.rhs)), 1e-30))
        velocity = self.fixed_velocity.copy()
        velocity[self.free] = state[:len(self.free)]
        pressure = np.r_[state[len(self.free):], 0.0]
        continuity = float(np.max(np.abs(self.divergence@velocity))/self.inflow)
        if (not np.isfinite(state).all() or linear > self.problem.linear_residual_limit or
                continuity > self.problem.linear_residual_limit):
            raise ArithmeticError(f"Stokes solve failed: linear={linear:g}, continuity={continuity:g}")
        dissipation = float(velocity@(stiffness@velocity))
        pressure_power = float(self.pressure_weights@pressure)
        rhs_adjoint = np.zeros_like(state)
        if objective == "dissipation":
            value = dissipation
            rhs_adjoint[:len(self.free)] = 2*(stiffness@velocity)[self.free]
            direct = np.asarray(self.drag_weights.T@(velocity**2)).ravel()
        else:
            value = pressure_power
            rhs_adjoint[len(self.free):] = self.pressure_weights[:-1]
            direct = np.zeros(self.cells)
        adjoint = factor.solve(rhs_adjoint, trans="T")
        adjoint_error = float(np.max(np.abs(matrix.T@adjoint-rhs_adjoint))/
                              max(np.max(np.abs(rhs_adjoint)), 1e-30))
        if not np.isfinite(adjoint).all() or adjoint_error > self.problem.linear_residual_limit:
            raise ArithmeticError(f"discrete adjoint failed: residual={adjoint_error:g}")
        implicit = np.asarray(self.free_weights.T@(adjoint[:len(self.free)]*velocity[self.free])).ravel()
        gradient = self.chain(projection_derivative, (direct-implicit)*dalpha)
        volume_gradient = self.chain(projection_derivative, np.full(self.cells, 1/self.cells))
        uc, vc = self.cell_velocity(velocity)
        speed = np.hypot(uc, vc)
        solid = rho < 0.1
        leakage = float(np.max(speed[solid])/max(np.max(speed), 1e-30)) if np.any(solid) else 0.0
        return Evaluation(value, gradient, float(np.mean(rho)), volume_gradient, rho,
                          velocity, pressure, linear, continuity, adjoint_error,
                          dissipation, pressure_power, float(np.mean(4*rho*(1-rho))), leakage)

    def cell_velocity(self, velocity):
        u = velocity[:self.nu_faces].reshape(self.ny, self.nx+1)
        v = velocity[self.nu_faces:].reshape(self.ny+1, self.nx)
        return ((u[:, :-1]+u[:, 1:])/2).ravel(), ((v[:-1]+v[1:])/2).ravel()

    def oc_candidate(self, x, evaluation, beta, move=0.15):
        """Moving-bound optimality-criteria step with physical-volume bisection.

        The caller evaluates and backtracks this proposal; it is never assumed
        to improve the objective. Positive sensitivities may move toward solid.
        """
        if not 0 < move <= 1:
            raise ValueError("move must be in (0,1]")
        x = self.enforce_passive(x)
        lower, upper = np.maximum(0, x-move), np.minimum(1, x+move)
        def candidate(multiplier):
            ratio = np.maximum(0, -evaluation.gradient)/(multiplier*
                    np.maximum(evaluation.volume_gradient, 1e-30))
            y = np.clip(np.maximum(x, 1e-3)*np.sqrt(ratio), lower, upper)
            return self.enforce_passive(y)
        low, high = 0.0, 1.0
        target = self.problem.volume_fraction
        for _ in range(100):
            if self.volume(candidate(high), beta)[0] <= target+1e-12:
                break
            high *= 2
        else:
            raise ArithmeticError("OC volume constraint could not be bracketed")
        for _ in range(55):
            mid = (low+high)/2
            if self.volume(candidate(mid), beta)[0] > target:
                low = mid
            else:
                high = mid
        return candidate(high)
