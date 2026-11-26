"""Polynomial system description and a parser for the repo's problems/*.m files."""
import random
import re
from dataclasses import dataclass, field

import sympy

from . import ff


@dataclass
class Problem:
    """A square or overdetermined system f_k(x; c) = 0.

    vars   : unknowns a1..an (sympy symbols)
    params : data coefficients c1..cm (sympy symbols)
    eqs    : sympy expressions, polynomial in vars with coefficients polynomial in params
    cfg    : generator options (hidden_var is 1-based, -1 = try all)
    """
    name: str
    vars: list
    params: list
    eqs: list
    cfg: dict = field(default_factory=dict)

    def __post_init__(self):
        self.n = len(self.vars)
        self.polys = []  # per equation: {exponent tuple: coefficient expr in params}
        for e in self.eqs:
            P = sympy.Poly(sympy.expand(e), *self.vars)
            self.polys.append({m: c for m, c in P.terms() if c != 0})
        self.supports = [sorted(p) for p in self.polys]
        # unique coefficient expressions, shared by all template entries
        self.coeffs, self._cid = [], {}
        self.coeff_ids = []
        for p in self.polys:
            ids = {}
            for m, c in p.items():
                key = sympy.srepr(c)
                if key not in self._cid:
                    self._cid[key] = len(self.coeffs)
                    self.coeffs.append(c)
                ids[m] = self._cid[key]
            self.coeff_ids.append(ids)

    def instance_mod_p(self, seed=0, p=ff.P):
        """Coefficient values mod p at a seeded random point of the data space."""
        rng = random.Random(seed)
        vals = {c: rng.randrange(1, p) for c in self.params}
        out = []
        for c in self.coeffs:
            P = sympy.Poly(c, *self.params) if self.params else sympy.Poly(c, sympy.Dummy())
            acc = 0
            for mon, q in P.terms():
                q = sympy.Rational(q)
                t = int(q.p) * ff.inv(int(q.q), p) % p
                for s, e in zip(P.gens, mon):
                    if e:
                        t = t * pow(vals.get(s, 0), e, p) % p
                acc = (acc + t) % p
            out.append(acc)
        return out


# ---------------------------------------------------------------------------
# MATLAB problem file parser
# ---------------------------------------------------------------------------

def _matlab_value(txt):
    txt = txt.strip()
    if txt.startswith('['):
        nums = re.findall(r'-?\d+', txt)
        return [int(x) for x in nums]
    try:
        return int(txt)
    except ValueError:
        return txt


def load_matlab_problem(path):
    """Parse problems/<name>.m (retrieve_eqs + retrieve_solver_cfg) into a Problem."""
    with open(path) as fh:
        src = fh.read()
    # join MATLAB line continuations
    src = re.sub(r'\.\.\.[^\n]*\n', ' ', src)
    m = re.search(r'function\s+eqs\s*=\s*retrieve_eqs\s*\(([^)]*)\)', src)
    if not m:
        raise ValueError(f'{path}: no retrieve_eqs function found')
    args = [a.strip() for a in m.group(1).split(',')]
    syms = {a: sympy.Symbol(a) for a in args}
    vars_ = [syms[a] for a in args if a.startswith('a')]
    params = [syms[a] for a in args if a.startswith('c')]
    eqs = {}
    for em in re.finditer(r'^\s*eqs\((\d+)\)\s*=\s*([^;]*);', src, re.M):
        expr = em.group(2).replace('.^', '^').replace('.*', '*').replace('./', '/').replace('^', '**')
        eqs[int(em.group(1))] = sympy.parse_expr(expr, local_dict=syms)
    cfg_raw = {k.lower(): _matlab_value(v) for k, v in re.findall(r'^\s*cfg\.(\w+)\s*=\s*(\[[^\]]*\]|[^;\n]*)', src, re.M)}
    cfg = {
        'hidden_var': cfg_raw.get('hiddenvarnum', -1),
        'poly_comb': cfg_raw.get('polycomb') or [],
        'size_of_combs': cfg_raw.get('sizeofcombs') or [],
        'rows_to_reduce': cfg_raw.get('noofrowstoreduce', 0),
    }
    name = re.sub(r'\.m$', '', path.split('/')[-1])
    return Problem(name, vars_, params, [eqs[k] for k in sorted(eqs)], cfg)
