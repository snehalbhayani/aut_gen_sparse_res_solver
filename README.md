# Automatic generator for a sparse resultant based polynomial solver

## Generating a solver (OFFLINE stage)
### Input
- A 'problem_name.m' file which returns a structure of solver configuration parameters and a function that returns a set of input polynomial equations
- The configuration parameters are housed in a matlab struct 'cfg'.
- The function that returns the input polynomials has the signature
> function eqs = retrieve_eqs(a1,a2,..,c1,c2,...)
- The configuration struct has the following fields
    - numOfCoeff
    - numOfVars
    - hiddenVarNum
    - sizeOfCombs or polyComb
    - noOfRowsToReduce
    - heurisiticTemplatesize
### Output 
- A 'solver.m' file in 'solvers/problem_name'
- Two other files are generated in 'solvers/problem_name' which are to be checked for debugging purposes.
    - eqs.txt
    - A MAPLE script which was executed for generating the solver
        
### Usage
- Execute
  > build_test_solver(`p1`, `p2`, `p3`)
  - `p1` is 1 if we want to generate a solver , 0 if we do not want to generate a solver
  - `p2` is 1 if we want to test a solver , 0 if we do not want to test a solver
  - `p3` is the number of random datapoints to be used for testing a solver
