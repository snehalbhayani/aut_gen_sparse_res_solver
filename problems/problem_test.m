function[eqsHandler, cfg] = problem_test()
%% Configuring the solver
cfg = retrieve_solver_cfg();
eqsHandler = @problem_test_eqs;
end
%% The polynomial system

%% The polynomial system
function cfg = retrieve_solver_cfg()
cfg = {};
cfg.numOfCoeff = 596;
cfg.numOfVars = 3;
% The index of the selected variable, x[i] or x_i in the extra polynomial,
% x_i - l.
% If set to -1, all variables will be tested one
% by one
cfg.hiddenVarNum = 1;
% (1) Either the size of polynomial combinations to be tested.
cfg.sizeOfCombs = [2];
% (2) Or the specific polynomial combination to be tested.
cfg.polyComb=[];
% The number of rows to be GJ eliminated to obtained a reduced input
% polynomial system as an input to the generator.
cfg.noOfRowsToReduce = 0;
% The heuristic size of the template. There is no theory to the best of our
% knowledge, governing the smallest template that can be generated.
% Hence one can start with a larger size and try to test by reducing the size
% of the template
cfg.heurisiticTemplatesize = 200;

cfg.eqsInst = @calib_3vars_exp;
end