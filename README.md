# Automatic generator for a sparse resultant based polynomial solver
## Required software
- Software: MATLAB and Maple.
- Currently supported MATLAB version: R2018a+
- Currently supported Maple version: 2018+
## Setup
- MATLAB installation should include the `symbolic maple toolbox` 
- The Maple installation has to be setup as the symbolic engine in MATLAB insallation
    - Older MATLAB versions connect to Maple through their `symbolic math toolboz`. 
    - But MATLAB 2018+ have a separate Maple toolbox to be setup.
- One quick way to check is by executing one of the following in the command MATLAB window
    > maple 
    - If Maple is connected, the command should open a GUI Maple interface.
- For more help on installing Maple toolbox for MATLAB, one can refer to <https://www.maplesoft.com/support/install/mtm11Install.html>.
   
## Generating a solver (OFFLINE stage)
### Input
- A `problem_name.m` file which returns a structure of solver configuration parameters and a function that returns a set of input polynomial equations
- The `problem_name.m` file is to be stored in `problems/` folder.
- The configuration parameters are housed in a matlab struct `cfg`.
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
- A `solver.m` file in 'solvers/problem_name'
- Two other files are generated in `solvers/problem_name` which are to be used for debugging purposes:
    - eqs.txt
    - A MAPLE script which was executed for generating the solver for problem `problem_name`
        
### Usage
- Navigate to the main folder of the generator
- Execute
  > build_test_solver(`p1`, `p2`, `p3`)
  - `p1` is 1 if we want to generate a solver , 0 if we do not want to generate a solver
  - `p2` is 1 if we want to test a solver , 0 if we do not want to test a solver
  - `p3` is the number of random datapoints to be used for testing a solver
- When prompted for the problem name, enter the value of `problem_name`

## Executing the solver (ONLINE stage)
- The solver for a problem problem_name is housed in `solvers/problem_name`.
- Execute 
    > build_test_solver(0, 1, `p`)
    - `p` is th number of random instances to be used for testing the solver
- When prompted for the problem name, enter the value of `problem_name`

## Further questions or comments
- Please write to snehal.bhayani@oulu.fi or snehalbhayani04@gmail.com

## Reference
- If you are using this generator software please cite the following:
<a id="1">[1]</a> 
Bhayani, S., Kukelova, Z., & Heikkilä, J. (2019). 
A sparse resultant based method for efficient minimal solvers. 
ArXiv, abs/1912.10268.
[PDF](https://arxiv.org/pdf/1912.10268.pdf)
