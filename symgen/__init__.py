"""symgen: deterministic, exact generator of sparse resultant minimal solvers.

    from symgen.problem import load_matlab_problem
    from symgen.template import search
    from symgen.codegen import matlab_solver

    pb = load_matlab_problem('problems/problem_relpose_6p_focal.m')
    T = search(pb)                       # exact template search + reduction
    open('solver.m', 'w').write(matlab_solver(pb, T))
"""
