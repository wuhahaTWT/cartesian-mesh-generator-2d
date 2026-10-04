"""Steady incompressible MAC Navier--Stokes--Brinkman and discrete adjoint.

Lengths use the design coordinates, inlet peak speed is 1, and mu is the
dimensionless viscous coefficient (default 1). With U_mean=2/3 and inlet width
w, gamma=Re*mu/(U_mean*w): gamma div(u tensor u)-mu Laplacian(u)+grad(p)+alpha u=0.
Thus Re=0 is an exact Stokes operator, not a small but finite Reynolds number.

Convection is the conservative centred MAC dual-volume flux: arithmetic
averages locate both transported and transporting velocities on a dual face.
This is the second-order Harlow--Welch stencil (1965, doi:10.1063/1.1761178).
Walls have zero normal flux. No artificial viscosity/limiter is added. Newton
convergence does NOT establish resolution or boundedness at high cell Re.

Each quadratic monomial is differentiated analytically. The factor passed to
the inherited adjoint is the FULL Jacobian at the accepted nonlinear root,
never an Oseen/Picard matrix. Filter/projection/material chain rules remain in
brinkman.py. Newton uses residual backtracking and adaptive continuation from
Stokes; failed attempts never replace the last converged continuation state.
"""
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu

from brinkman import StokesBrinkman


class NavierStokesBrinkman(StokesBrinkman):
    def __init__(self, problem):
        super().__init__(problem)
        self.inertia = problem.reynolds*problem.viscosity/((2/3)*self.aperture_width)
        self.newton_history = []
        self._convection_stencil()

    def _convection_stencil(self):
        rows, aa, bb, cc = [], [], [], []
        def product(row, first, second, weight):
            for a in first:
                for b in second:
                    rows.append(row); aa.append(a); bb.append(b)
                    cc.append(weight/(len(first)*len(second)))
        for j in range(self.ny):
            for i in range(1, self.nx):
                row = self.u(i, j)
                for delta, sign in ((0, 1), (-1, -1)):
                    pair = [self.u(i+delta, j), self.u(i+delta+1, j)]
                    product(row, pair, pair, sign*self.dy)
                if j+1 < self.ny:
                    product(row, [self.v(i-1, j+1), self.v(i, j+1)],
                            [self.u(i, j), self.u(i, j+1)], self.dx)
                if j > 0:
                    product(row, [self.v(i-1, j), self.v(i, j)],
                            [self.u(i, j-1), self.u(i, j)], -self.dx)
        for j in range(1, self.ny):
            for i in range(self.nx):
                row = self.v(i, j)
                for delta, sign in ((0, 1), (-1, -1)):
                    pair = [self.v(i, j+delta), self.v(i, j+delta+1)]
                    product(row, pair, pair, sign*self.dx)
                if i+1 < self.nx:
                    product(row, [self.u(i+1, j-1), self.u(i+1, j)],
                            [self.v(i, j), self.v(i+1, j)], self.dy)
                if i > 0:
                    product(row, [self.u(i, j-1), self.u(i, j)],
                            [self.v(i-1, j), self.v(i, j)], -self.dy)
        self._cr, self._ca, self._cb = (np.asarray(a, dtype=int) for a in (rows, aa, bb))
        self._cc = np.asarray(cc)

    def convection(self, velocity, jacobian=True):
        r, a, b, c = self._cr, self._ca, self._cb, self._cc
        value = np.bincount(r, weights=c*velocity[a]*velocity[b], minlength=self.velocities)
        if not jacobian:
            return value
        derivative = sparse.coo_matrix((np.r_[c*velocity[b], c*velocity[a]],
                    (np.r_[r, r], np.r_[a, b])), shape=(self.velocities,)*2).tocsc()
        return value, derivative

    def _solve(self, stiffness):
        state, stokes, factor, residual = super()._solve(stiffness)
        self.newton_history = [dict(inertia=0., residual=residual, iterations=0, accepted=True)]
        if self.inertia == 0:
            return state, stokes, factor, residual
        n = len(self.free)
        # Residual: infinity norm of [integrated momentum; -integrated mass]
        # divided by max(||eliminated Dirichlet rhs||_inf,1e-30), exactly the
        # legacy algebraic normalization. Default <=1e-11, plus inherited
        # max|D u|/Q_in <=1e-8. These are equation, NOT physical-error, gates.
        scale = max(np.max(np.abs(self.rhs)), 1e-30)
        def equation(z, coefficient, with_jacobian=True):
            velocity = self.fixed_velocity.copy()
            velocity[self.free] = z[:n]
            if with_jacobian:
                flux, jac = self.convection(velocity)
                block = sparse.bmat([[jac[self.free, :][:, self.free],
                         sparse.csc_matrix((n, self.cells-1))],
                         [sparse.csc_matrix((self.cells-1, n)),
                          sparse.csc_matrix((self.cells-1, self.cells-1))]], format="csc")
            else:
                flux = self.convection(velocity, False)
            f = stokes@z-self.rhs
            f[:n] += coefficient*flux[self.free]
            return (f, stokes+coefficient*block) if with_jacobian else f
        coefficient, step = 0., self.inertia
        attempts = 0
        while coefficient < self.inertia:
            target = min(self.inertia, coefficient+step)
            trial = state.copy()
            accepted = False
            for iteration in range(self.problem.newton_iterations+1):
                f, matrix = equation(trial, target)
                norm = float(np.max(np.abs(f))/scale)
                if np.isfinite(norm) and norm <= self.problem.newton_residual_limit:
                    accepted = True
                    break
                if not np.isfinite(norm) or iteration == self.problem.newton_iterations:
                    break
                delta = splu(matrix).solve(-f)
                damping = 1.
                for _ in range(20):
                    proposal = trial+damping*delta
                    error = np.max(np.abs(equation(proposal, target, False)))/scale
                    if np.isfinite(error) and error <= (1-1e-4*damping)*norm:
                        trial = proposal
                        break
                    damping *= .5
                else:
                    break
            self.newton_history.append(dict(inertia=target, residual=norm,
                                            iterations=iteration, accepted=accepted))
            attempts += 1
            if accepted:
                state, coefficient = trial, target
                step *= 1.5
            else:
                step *= .5
                if step < self.inertia*1e-6 or attempts > 80:
                    raise ArithmeticError(f"Newton continuation failed at inertia={target:g}, residual={norm:g}")
        # Refactor at the actual converged root for the transpose adjoint.
        return state, matrix, splu(matrix), norm
