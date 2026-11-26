"""Floating point evaluation of a template (reference solver + accuracy check)."""
import numpy as np
import sympy


class SolverLayout:
    """Index bookkeeping shared by the NumPy solver and the MATLAB code generator."""

    def __init__(self, problem, T):
        self.pb, self.T = problem, T
        n, h = problem.n, T.h
        eh = tuple(1 if j == h else 0 for j in range(n))
        self.E, self.G = list(T.E), list(T.G)
        Ei = {e: i for i, e in enumerate(self.E)}
        Gi = {g: i for i, g in enumerate(self.G)}
        self.nE, self.nG = len(self.E), len(self.G)
        # A2 (nG x nG) and A1 (nG x nE) entries: (row, col, coefficient id)
        self.A2, self.A1 = [], []
        for i, (k, t) in enumerate(T.top):
            for a, cid in problem.coeff_ids[k].items():
                m = tuple(x + y for x, y in zip(t, a))
                if m in Gi:
                    self.A2.append((i, Gi[m], cid))
                else:
                    self.A1.append((i, Ei[m], cid))
        # action matrix rows: lambda * v_e = v_{e + e_h}
        self.X_unit, self.X_from_K = [], []   # (row in X, col in E) / (row in X, row in G)
        for i, e in enumerate(self.E):
            m = tuple(x + y for x, y in zip(e, eh))
            if m in Ei:
                self.X_unit.append((i, Ei[m]))
            else:
                self.X_from_K.append((i, Gi[m]))
        # recovery: monomial index in W = [v_E; v_G]
        widx = {**Ei, **{g: self.nE + i for g, i in Gi.items()}}
        self.recovery = {i: [(widx[e], widx[xe]) for e, xe in pairs] for i, pairs in T.recovery.items()}
        self.used_coeffs = sorted({c for _, _, c in self.A1 + self.A2})


class NumericSolver:
    def __init__(self, problem, T):
        self.pb, self.L = problem, SolverLayout(problem, T)
        self.h = T.h
        self._coef = sympy.lambdify(problem.params, [problem.coeffs[c] for c in self.L.used_coeffs], 'numpy')
        self._cpos = {c: i for i, c in enumerate(self.L.used_coeffs)}
        self._eqs = sympy.lambdify(problem.vars + problem.params, problem.eqs, 'numpy')

    def solve(self, data):
        L = self.L
        k = np.array(self._coef(*data), dtype=complex)
        A2 = np.zeros((L.nG, L.nG), complex)
        A1 = np.zeros((L.nG, L.nE), complex)
        for i, j, c in L.A2:
            A2[i, j] = k[self._cpos[c]]
        for i, j, c in L.A1:
            A1[i, j] = k[self._cpos[c]]
        K = np.linalg.solve(A2, A1)
        X = np.zeros((L.nE, L.nE), complex)
        for i, j in L.X_unit:
            X[i, j] = 1
        for i, g in L.X_from_K:
            X[i] = -K[g]
        lam, V = np.linalg.eig(X)
        ok = np.isfinite(lam)
        lam, V = lam[ok], V[:, ok]
        W = np.vstack([V, -K @ V])
        sols = np.zeros((len(lam), self.pb.n), complex)
        sols[:, self.h] = lam
        for i, pairs in L.recovery.items():
            P = np.array([p for p, _ in pairs])
            Q = np.array([q for _, q in pairs])
            b = np.argmax(np.abs(W[P]), axis=0)
            cols = np.arange(W.shape[1])
            sols[:, i] = W[Q[b], cols] / W[P[b], cols]
        return sols

    def residuals(self, sols, data):
        out = []
        for x in sols:
            r = np.abs(np.array(self._eqs(*x, *data), dtype=complex))
            out.append(np.max(r) / max(np.linalg.norm(x), 1e-300))
        return np.array(out)


def evaluate(problem, T, trials=50, seed=1, tol=1e-6):
    """Accuracy statistics on random (Gaussian) data."""
    rng = np.random.default_rng(seed)
    S = NumericSolver(problem, T)
    allres, counts = [], []
    for _ in range(trials):
        data = rng.standard_normal(len(problem.params))
        sols = S.solve(data)
        res = S.residuals(sols, data)
        allres.append(res)
        counts.append(int(np.sum(res < tol)))
    allres = np.concatenate(allres)
    return {
        'median_log10_residual': float(np.median(np.log10(allres + 1e-300))),
        'mean_log10_residual': float(np.mean(np.log10(allres + 1e-300))),
        'good_solutions_per_instance': float(np.mean(counts)),
        'eigen_size': S.L.nE,
    }
