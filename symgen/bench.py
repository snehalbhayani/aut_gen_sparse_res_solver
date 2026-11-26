"""Run a MATLAB solver.m in Octave on random data and measure accuracy in Python.

Works for both the Maple-generated solvers in solvers/ and the symgen ones,
since both return rows [a_h, other unknowns..., lambda].
"""
import os
import subprocess
import tempfile

import numpy as np

from .numeric import NumericSolver  # noqa: F401  (re-export for convenience)
import sympy


def run_octave_solver(folder, data, timeout=1200):
    with tempfile.TemporaryDirectory() as tmp:
        dfile = os.path.join(tmp, 'data.txt')
        ofile = os.path.join(tmp, 'out.txt')
        np.savetxt(dfile, data, fmt='%.17g')
        # Octave's rref errors on 0-row input (MATLAB returns empty); the Maple
        # solvers call rref(M(1:0,:)).  Shim it for benchmarking.
        with open(os.path.join(tmp, 'rref.m'), 'w') as fh:
            fh.write("function R = rref(A)\nif isempty(A), R = A; return; end\n"
                     "error('rref shim: only empty input supported');\nend\n")
        script = f"""
warning('off', 'all');
addpath('{os.path.abspath(folder)}');
addpath('{tmp}');
D = dlmread('{dfile}');
fid = fopen('{ofile}', 'w');
t = 0;
for i = 1:size(D,1)
  tic; S = solver(D(i,:)'); t = t + toc;
  for j = 1:size(S,1)
    fprintf(fid, '%d', i);
    fprintf(fid, ' %.17g %.17g', [real(S(j,:)); imag(S(j,:))]);
    fprintf(fid, '\\n');
  end
end
fclose(fid);
printf('%.6g\\n', t / size(D,1));
"""
        r = subprocess.run(['octave', '--no-gui', '--quiet', '--eval', script],
                           capture_output=True, text=True, timeout=timeout)
        if r.returncode != 0:
            raise RuntimeError(r.stderr[-2000:])
        sec = float(r.stdout.strip().splitlines()[-1])
        rows = np.loadtxt(ofile, ndmin=2) if os.path.getsize(ofile) else np.zeros((0, 1))
    sols = {}
    for row in rows:
        i = int(row[0]) - 1
        v = row[1::2] + 1j * row[2::2]
        sols.setdefault(i, []).append(v)
    return sols, sec


def octave_accuracy(problem, folder, hidden_var, trials=30, seed=1, tol=1e-6):
    rng = np.random.default_rng(seed)
    data = rng.standard_normal((trials, len(problem.params)))
    sols, sec = run_octave_solver(folder, data)
    f = sympy.lambdify(problem.vars + problem.params, problem.eqs, 'numpy')
    h = hidden_var - 1
    order = [h] + [i for i in range(problem.n) if i != h]
    res, counts = [], []
    for i in range(trials):
        c = 0
        for v in sols.get(i, []):
            x = np.zeros(problem.n, complex)
            x[order] = v[:problem.n]
            if not np.all(np.isfinite(x)):
                continue
            r = np.max(np.abs(np.array(f(*x, *data[i]), dtype=complex))) / max(np.linalg.norm(x), 1e-300)
            res.append(r)
            c += r < tol
        counts.append(c)
    res = np.array(res)
    return {'median_log10_residual': float(np.median(np.log10(res + 1e-300))),
            'good_solutions_per_instance': float(np.mean(counts)),
            'seconds_per_solve_octave': sec}
