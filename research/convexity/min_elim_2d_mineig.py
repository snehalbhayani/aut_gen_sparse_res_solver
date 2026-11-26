"""Part 1b probe: smallest |G| when the eigen size is the minimum (9 = 8 real + 1 ghost)."""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'min_elim_2d.py')).read().split('best = None')[0])
best = None
for W in [(1, 1), (2, 1), (1, 2), (1, 3), (3, 1), (2, 3), (3, 2)]:
    tb = TemplateBuilder(pb, 0, weights=W)
    for kx, ky in itertools.product(range(3, 9), range(3, 10)):
        B = [(a, b) for a in range(kx + 1) for b in range(ky + 1)]
        rep = tb.repair(*tb.from_basis(B))
        if rep is None:
            continue
        top, E = tb.reduce(*rep)
        cols = tb.columns(top, E); G = sorted(cols - set(E))
        if len(E) == 9 and (best is None or len(G) < best[0]):
            best = (len(G), (kx, ky), W, sorted(cols))
print(f'eigen size 9: smallest |G| found = {best[0]} (start box {best[1]}, order {best[2]}); |C| = {len(best[3])}')
ys = range(max(p[1] for p in best[3]), -1, -1)
for y in ys:
    print('   ' + ''.join('#' if (x, y) in best[3] else '.' for x in range(max(p[0] for p in best[3]) + 1)))
