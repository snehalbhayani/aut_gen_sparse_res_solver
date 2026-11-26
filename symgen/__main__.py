"""symgen command line: generate a sparse resultant solver from problems/<name>.m

    python3 -m symgen problems/problem_relpose_6p_focal.m [options]

Writes solvers_symgen/<name>/solver.m (same calling convention as the Maple
generated solvers) and template.json (the template, for inspection).
"""
import argparse
import json
import os
import time

from .codegen import matlab_solver, template_json
from .numeric import evaluate
from .problem import load_matlab_problem
from .template import search


def _ints(s):
    return [int(x) for x in s.replace(';', ',').split(',') if x.strip()]


def main(argv=None):
    ap = argparse.ArgumentParser(prog='symgen', description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('problem', help='path to problems/<name>.m, or module:function returning a Problem')
    ap.add_argument('--out', default='solvers_symgen', help='output root folder')
    ap.add_argument('--hidden', type=int, help='hidden variable (1-based, -1 = try all); default from cfg')
    ap.add_argument('--comb', type=_ints, help='polynomial combination, e.g. 1,11 (m+1 = extra polynomial)')
    ap.add_argument('--sizes', type=_ints, help='try all combinations of these sizes, e.g. 2,3')
    ap.add_argument('--extra', choices=['segment', 'simplex', 'both'], default='both',
                    help='Newton polytope used for the extra polynomial x_h - lambda')
    ap.add_argument('--eigen-target', type=int, help='stop eigen basis reduction at this size')
    ap.add_argument('--fastest', action='store_true',
                    help='sweep eigen size vs elimination size, time each solver, keep the fastest')
    ap.add_argument('--keep', type=int, default=6, help='number of candidate templates to fully reduce')
    ap.add_argument('--seed', type=int, default=0, help='seed of the GF(p) instance used for exact rank tests')
    ap.add_argument('--test', type=int, default=30, help='random instances for the accuracy check (0 = skip)')
    ap.add_argument('-v', '--verbose', action='store_true', help='log every candidate template')
    ap.add_argument('--no-cse', action='store_true', help='do not extract common subexpressions')
    a = ap.parse_args(argv)

    if a.problem.endswith('.m'):
        pb = load_matlab_problem(a.problem)
    else:  # python problem: package.module:function returning a Problem
        import importlib
        mod, fn = a.problem.split(':')
        pb = getattr(importlib.import_module(mod), fn)()
    cfg = {'extra': a.extra}
    if a.hidden is not None:
        cfg['hidden_var'] = a.hidden
    if a.comb:
        cfg['poly_comb'] = a.comb
    if a.sizes:
        cfg['poly_comb'], cfg['size_of_combs'] = [], a.sizes
    print(f'{pb.name}: {pb.n} unknowns, {len(pb.eqs)} equations, {len(pb.params)} data coefficients')
    t0 = time.time()
    if a.eigen_target:
        cfg['eigen_target'] = a.eigen_target
    if a.fastest:
        from .fastest import sweep
        best, _ = sweep(pb, cfg, seed=a.seed)
        T = best['template']
    else:
        T = search(pb, cfg, seed=a.seed, keep=a.keep, verbose=a.verbose)
    print(f'template: {T.size[1]} x {T.size[1]} elimination + {T.size[0]} x {T.size[0]} eigenvalue problem '
          f'(hidden a{T.h + 1}), found in {time.time() - t0:.1f}s')
    folder = os.path.join(a.out, pb.name)
    os.makedirs(folder, exist_ok=True)
    with open(os.path.join(folder, 'solver.m'), 'w') as fh:
        fh.write(matlab_solver(pb, T, cse=not a.no_cse))
    with open(os.path.join(folder, 'template.json'), 'w') as fh:
        fh.write(template_json(pb, T))
    print(f'wrote {folder}/solver.m and template.json')
    if a.test:
        stats = evaluate(pb, T, trials=a.test)
        print('accuracy (NumPy reference):', json.dumps(stats))


if __name__ == '__main__':
    main()
