"""Randomized template reduction: how far is symgen's deterministic greedy from optimal?

symgen removes columns, rows and eigen monomials in one fixed priority order.  Here the
priority is a random permutation per run (all reductions stay exact over GF(p)), started
from every valid Minkowski-sum basis.  Reports the best |G| at the minimal eigen size.

    python3 research/convexity/random_reduce.py problem_9pt2radial 3 27 200 [--no-recovery]

--no-recovery drops the variable-recovery requirement (pairs e, e + e_i in the template),
which is the setting of the lower bounds in research/proof/lower_bounds.md.
"""
import os, random, sys, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from symgen.problem import load_matlab_problem
from symgen.template import TemplateBuilder, candidate_bases, default_combs


class RandomBuilder(TemplateBuilder):
    def __init__(self, *a, rng=None, recovery=True, **k):
        super().__init__(*a, **k)
        self.key = {}
        self.rng = rng
        self.need_recovery = recovery

    def _prio(self, m):
        if m not in self.key:
            self.key[m] = self.rng.random()
        return self.key[m]

    def recoverable(self, E, cols):
        return super().recoverable(E, cols) if self.need_recovery else True


def main():
    name, h, N, runs = sys.argv[1], int(sys.argv[2]) - 1, int(sys.argv[3]), int(sys.argv[4])
    recovery = '--no-recovery' not in sys.argv
    pb = load_matlab_problem(f'problems/{name}.m')
    cfg = dict(pb.cfg)
    cfg['size_of_combs'] = cfg.get('size_of_combs') or [pb.n, pb.n + 1]
    base = RandomBuilder(pb, h, rng=random.Random(0), recovery=recovery)
    starts = []
    for desc, B in candidate_bases(pb, h, default_combs(pb, cfg)):
        top, E = base.from_basis(B)
        rep = base.repair(top, E)
        if rep is not None:
            starts.append((desc, rep))
    print(f'{name}, hidden a{h+1}, recovery {"on" if recovery else "off"}: {len(starts)} valid start bases', flush=True)
    best, t0 = None, time.time()
    rng = random.Random(1)
    for run in range(runs):
        desc, (top, E) = starts[run % len(starts)]
        tb = RandomBuilder(pb, h, rng=random.Random(rng.random()), recovery=recovery)
        top2, E2 = tb.reduce(top, E)
        T = tb.finalize(top2, E2)
        if len(T.E) == N and (best is None or len(T.G) < best[0]):
            best = (len(T.G), desc, run)
            print(f'  run {run:4d} ({time.time()-t0:6.0f}s): eigen {len(T.E)}, |G| = {len(T.G)}  from {desc}', flush=True)
    print(f'best |G| at eigen {N}: {best[0] if best else None} over {runs} runs')


if __name__ == '__main__':
    main()
