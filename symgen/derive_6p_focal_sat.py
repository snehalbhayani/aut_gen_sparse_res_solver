"""Derive a 6p-focal formulation without the 3 spurious solutions at a3 = infinity.

Original system (problems/problem_relpose_6p_focal.m):
    F = M(c1..9) + a1*M(c10..18) + a2*M(c19..27)     (row major 3x3 blocks)
    Q = diag(1, 1, a3)
    det F = 0,   E = 2 F Q F' Q F - trace(F Q F' Q) F = 0   (9 equations)

E = W0 + a3 W1 + a3^2 W2 with W2 = F33 (2 f3 r3 - F33 F).  At a3 = infinity
only W2 has to vanish, and W2 = 0 together with det F = 0 has exactly 3
solutions (F33 = 0 is a line, det F = 0 a cubic).  Every template built from
the original equations therefore counts 15 + 3 = 18 solutions.

The leading coefficients satisfy the syzygies (i, j in {1, 2})
    F33 W2_i3 - F_i3 W2_33 = 0
    F33 W2_3j - F_3j W2_33 = 0
    F33 W2_ij - 2 F_3j W2_i3 + F_ij W2_33 = 0
so the same combinations of the E_ij have a3-degree <= 1.  They vanish at every
real solution, but not at the points at infinity: this is a symbolic
saturation of the ideal w.r.t. the hyperplane at infinity in a3, with
coefficients that are explicit polynomials in the data.

Run:  python3 -m symgen.derive_6p_focal_sat   -> problems/problem_relpose_6p_focal_sat.m
      python3 -m symgen symgen.derive_6p_focal_sat:make_problem --comb 10,11,19
The .m file has every coefficient expanded in the data (for the MATLAB
pipeline).  make_problem() keeps the structure instead: the coefficients of
E_ij are intermediate symbols and the new equations are short bilinear
expressions in them, which makes the generated solver much cheaper.
"""
import os

import sympy
from sympy.printing.octave import octave_code

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def derive():
    a1, a2, a3 = sympy.symbols('a1:4')
    c = sympy.symbols('c1:28')
    M = lambda k: sympy.Matrix(3, 3, list(c[9 * k: 9 * k + 9]))
    F = M(0) + a1 * M(1) + a2 * M(2)
    Q = sympy.diag(1, 1, a3)
    E = (2 * F * Q * F.T * Q * F - (F * Q * F.T * Q).trace() * F).applyfunc(sympy.expand)
    det = sympy.expand(F.det())
    W2 = E.applyfunc(lambda e: sympy.expand(e).coeff(a3, 2))
    # sanity check of the leading form
    assert all(sympy.expand(W2[i, j] - F[2, 2] * (2 * F[i, 2] * F[2, j] - F[2, 2] * F[i, j])) == 0
               for i in range(3) for j in range(3))
    new = []
    for i in range(2):
        new.append(F[2, 2] * E[i, 2] - F[i, 2] * E[2, 2])
        new.append(F[2, 2] * E[2, i] - F[2, i] * E[2, 2])
    for i in range(2):
        for j in range(2):
            new.append(F[2, 2] * E[i, j] - 2 * F[2, j] * E[i, 2] + F[i, j] * E[2, 2])
    new = [sympy.expand(g) for g in new]
    for g in new:
        assert sympy.degree(g, a3) <= 1, 'a3^2 term did not cancel'
    orig = [E[i, j] for i in range(3) for j in range(3)] + [det]
    return [a1, a2, a3], list(c), orig, new


def make_problem():
    """Structured Problem: E_ij coefficients as intermediate (aux) symbols."""
    from .problem import Problem
    a1, a2, a3 = sympy.symbols('a1:4')
    vars_ = [a1, a2, a3]
    c = sympy.symbols('c1:28')
    M = lambda k: sympy.Matrix(3, 3, list(c[9 * k: 9 * k + 9]))
    F = M(0) + a1 * M(1) + a2 * M(2)
    Q = sympy.diag(1, 1, a3)
    E = (2 * F * Q * F.T * Q * F - (F * Q * F.T * Q).trace() * F).applyfunc(sympy.expand)
    aux, Es = [], sympy.zeros(3, 3)
    for i in range(3):
        for j in range(3):
            P = sympy.Poly(E[i, j], *vars_)
            acc = 0
            for k, (mon, coef) in enumerate(P.terms()):
                s = sympy.Symbol(f'e{i + 1}{j + 1}_{k + 1}')
                aux.append((s, coef))
                acc += s * sympy.prod([v ** p for v, p in zip(vars_, mon)])
            Es[i, j] = acc
    new = []
    for i in range(2):
        new.append(F[2, 2] * Es[i, 2] - F[i, 2] * Es[2, 2])
        new.append(F[2, 2] * Es[2, i] - F[2, i] * Es[2, 2])
    for i in range(2):
        for j in range(2):
            new.append(F[2, 2] * Es[i, j] - 2 * F[2, j] * Es[i, 2] + F[i, j] * Es[2, 2])
    eqs = [Es[i, j] for i in range(3) for j in range(3)] + [sympy.expand(F.det())] + [sympy.expand(g) for g in new]
    cfg = {'hidden_var': 2, 'poly_comb': [10, 11, 19], 'size_of_combs': [], 'rows_to_reduce': 0}
    return Problem('problem_relpose_6p_focal_sat', vars_, list(c), eqs, cfg, aux)


def write_problem(path):
    vars_, params, orig, new = derive()
    eqs = orig + new
    args = ','.join(str(s) for s in vars_ + params)
    lines = ['function[eqsHandler, cfg] = problem_relpose_6p_focal_sat()',
             '% 6p relative pose with unknown focal length, with the 3 spurious solutions at',
             '% a3 = infinity removed by symbolic saturation (see symgen/derive_6p_focal_sat.py).',
             '% eqs 1-9: trace constraint, eq 10: det(F), eqs 11-18: a3-degree-1 combinations.',
             'cfg = retrieve_solver_cfg();', 'eqsHandler = @retrieve_eqs;', 'end',
             f'function eqs = retrieve_eqs({args})']
    for k, e in enumerate(eqs):
        lines.append(f'eqs({k + 1}) = {octave_code(e)};')
    lines += ['end', 'function cfg = retrieve_solver_cfg()', 'cfg = {};',
              f'cfg.numOfCoeff = {len(params)};', f'cfg.numOfVars = {len(vars_)};',
              'cfg.hiddenVarNum = 2;', 'cfg.sizeOfCombs = [];', 'cfg.polyComb = [11;19];',
              'cfg.noOfRowsToReduce = 0;', 'cfg.heurisiticTemplatesize = 100;', 'end']
    with open(path, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    out = os.path.join(ROOT, 'problems', 'problem_relpose_6p_focal_sat.m')
    write_problem(out)
    print('wrote', out)
