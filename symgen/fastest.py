"""Trade eigenvalue-problem size for elimination size and keep the fastest solver.

Each l-row dropped by the eigen basis reduction makes the eigenvalue problem
one smaller but adds one column to A2 (the matrix that is inverted).  Stopping
the reduction at an `eigen_target` gives the whole family of valid templates
between "largest eigenvalue problem, smallest A2" and "smallest eigenvalue
problem, largest A2".  This module sweeps that family, times every generated
solver (Octave if available, otherwise a flop model) and returns the fastest.

A larger eigenvalue problem than the number of solutions returns extra
(spurious) eigenpairs; they are reported in the table and can be discarded by
their residual.
"""
import os
import shutil
import tempfile

from .codegen import matlab_solver
from .numeric import evaluate
from .template import search


def flops(nE, nG):
    """Rough cost: LU solve of A2 with nE right hand sides + nonsymmetric eig (~10 n^3)."""
    return 2 / 3 * nG ** 3 + 2 * nG ** 2 * nE + 10 * nE ** 3


def sweep(problem, cfg=None, targets=None, steps=8, trials=200, seed=0, log=print):
    cfg = dict(cfg or {})
    base = search(problem, dict(cfg), seed=seed, verbose=False)
    nmin = base.size[0]
    # largest sensible eigen size: the unreduced l-rows of the best candidate
    nmax = base.info['initial'][0]
    if targets is None:
        step = max(1, (nmax - nmin) // steps)
        targets = sorted(set(list(range(nmin, nmax + 1, step)) + [nmax]))
    have_octave = shutil.which('octave') is not None
    rows = []
    for t in targets:
        T = base if t == nmin else search(problem, dict(cfg, eigen_target=t), seed=seed, verbose=False)
        row = {'eigen': T.size[0], 'elim': T.size[1], 'flops': flops(*T.size), 'template': T}
        acc = evaluate(problem, T, trials=30)
        row['good_solutions'] = acc['good_solutions_per_instance']
        row['median_log10_residual'] = acc['median_log10_residual']
        if have_octave:
            from .bench import octave_accuracy
            with tempfile.TemporaryDirectory() as d:
                with open(os.path.join(d, 'solver.m'), 'w') as fh:
                    fh.write(matlab_solver(problem, T))
                r = octave_accuracy(problem, d, T.h + 1, trials=trials)
            row['ms_octave'] = 1e3 * r['seconds_per_solve_octave']
        rows.append(row)
        log('eigen %3d  elim %3d  flops %8.0fk  good sols %5.2f  med log10 res %6.2f%s' % (
            row['eigen'], row['elim'], row['flops'] / 1e3, row['good_solutions'], row['median_log10_residual'],
            '  %.3f ms' % row['ms_octave'] if have_octave else ''))
    # prefer templates that still find (nearly) all solutions of the minimal one
    ok = [r for r in rows if r['good_solutions'] >= rows[0]['good_solutions'] - 0.5] or rows
    key = (lambda r: r['ms_octave']) if have_octave else (lambda r: r['flops'])
    return min(ok, key=key), rows
