# Automatic generator for a sparse resultant based polynomial solver

## Generating a solver
- A file containing input polynomials in symbolic form
- Variables are labeled as `a1`, `a2`, ... 
- Coefficients are labeled as `c1`, `c2`, `c3`, ... 
- > build_test_solver(`p1`, `p2`, `p3`)
- - `p1` is 1 if we want to generate a solver , 0 if we do not want to generate a solver
- - `p2` is 1 if we want to test a solver , 0 if we do not want to test a solver
- - `p3` is the number of random datapoints to be used for testing a solver
