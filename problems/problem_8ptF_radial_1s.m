function[eqsHandler, cfg] = problem_9pt2radial()
%% Configuring the solver
cfg = retrieve_solver_cfg();
% If the polynomial system is huge, this function can be
% extracted out in a separate file and the file can be loaded here to
% obtain the input polynomial system.s
eqsHandler = @retrieve_eqs;
end
%% The polynomial system of the solver
% The parameters are the variables and the coefficients. The variables have
% to be labelled as 'a1', 'a2', ... and the coefficients are labelled as
% 'c1', 'c2', ... 
function eqs = retrieve_eqs(a1,a2,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,c31,c32) 
eqs(1) = a1^3*a2^3*c9*c13*c25 - a1^3*a2^3*c9*c17*c21 - a1^3*a2^2*c1*c25 + a1^3*a2^2*c5*c21 + a1^3*a2^2*c9*c13*c27 + a1^3*a2^2*c9*c15*c25 - a1^3*a2^2*c9*c17*c23 - a1^3*a2^2*c9*c19*c21 + a1^3*a2^2*c11*c13*c25 - a1^3*a2^2*c11*c17*c21 - a1^3*a2*c1*c27 - a1^3*a2*c3*c25 + a1^3*a2*c5*c23 + a1^3*a2*c7*c21 + a1^3*a2*c9*c15*c27 - a1^3*a2*c9*c19*c23 + a1^3*a2*c11*c13*c27 + a1^3*a2*c11*c15*c25 - a1^3*a2*c11*c17*c23 - a1^3*a2*c11*c19*c21 - a1^3*c3*c27 + a1^3*c7*c23 + a1^3*c11*c15*c27 - a1^3*c11*c19*c23 + a1^2*a2^3*c9*c13*c26 + a1^2*a2^3*c9*c14*c25 - a1^2*a2^3*c9*c17*c22 - a1^2*a2^3*c9*c18*c21 + a1^2*a2^3*c10*c13*c25 - a1^2*a2^3*c10*c17*c21 + a1^2*a2^2*c1*c17 - a1^2*a2^2*c1*c26 - a1^2*a2^2*c2*c25 - a1^2*a2^2*c5*c13 + a1^2*a2^2*c5*c22 + a1^2*a2^2*c6*c21 + a1^2*a2^2*c9*c13*c28 + a1^2*a2^2*c9*c14*c27 + a1^2*a2^2*c9*c15*c26 + a1^2*a2^2*c9*c16*c25 - a1^2*a2^2*c9*c17*c24 - a1^2*a2^2*c9*c18*c23 - a1^2*a2^2*c9*c19*c22 - a1^2*a2^2*c9*c20*c21 + a1^2*a2^2*c10*c13*c27 + a1^2*a2^2*c10*c15*c25 - a1^2*a2^2*c10*c17*c23 - a1^2*a2^2*c10*c19*c21 + a1^2*a2^2*c11*c13*c26 + a1^2*a2^2*c11*c14*c25 - a1^2*a2^2*c11*c17*c22 - a1^2*a2^2*c11*c18*c21 + a1^2*a2^2*c12*c13*c25 - a1^2*a2^2*c12*c17*c21 + a1^2*a2*c1*c19 - a1^2*a2*c1*c28 - a1^2*a2*c2*c27 + a1^2*a2*c3*c17 - a1^2*a2*c3*c26 - a1^2*a2*c4*c25 - a1^2*a2*c5*c15 + a1^2*a2*c5*c24 + a1^2*a2*c6*c23 - a1^2*a2*c7*c13 + a1^2*a2*c7*c22 + a1^2*a2*c8*c21 + a1^2*a2*c9*c15*c28 + a1^2*a2*c9*c16*c27 - a1^2*a2*c9*c19*c24 - a1^2*a2*c9*c20*c23 + a1^2*a2*c10*c15*c27 - a1^2*a2*c10*c19*c23 + a1^2*a2*c11*c13*c28 + a1^2*a2*c11*c14*c27 + a1^2*a2*c11*c15*c26 + a1^2*a2*c11*c16*c25 - a1^2*a2*c11*c17*c24 - a1^2*a2*c11*c18*c23 - a1^2*a2*c11*c19*c22 - a1^2*a2*c11*c20*c21 + a1^2*a2*c12*c13*c27 + a1^2*a2*c12*c15*c25 - a1^2*a2*c12*c17*c23 - a1^2*a2*c12*c19*c21 + a1^2*c3*c19 - a1^2*c3*c28 - a1^2*c4*c27 - a1^2*c7*c15 + a1^2*c7*c24 + a1^2*c8*c23 + a1^2*c11*c15*c28 + a1^2*c11*c16*c27 - a1^2*c11*c19*c24 - a1^2*c11*c20*c23 + a1^2*c12*c15*c27 - a1^2*c12*c19*c23 + a1*a2^3*c9*c14*c26 - a1*a2^3*c9*c18*c22 + a1*a2^3*c10*c13*c26 + a1*a2^3*c10*c14*c25 - a1*a2^3*c10*c17*c22 - a1*a2^3*c10*c18*c21 + a1*a2^2*c1*c18 + a1*a2^2*c2*c17 - a1*a2^2*c2*c26 - a1*a2^2*c5*c14 - a1*a2^2*c6*c13 + a1*a2^2*c6*c22 + a1*a2^2*c9*c14*c28 + a1*a2^2*c9*c16*c26 - a1*a2^2*c9*c18*c24 - a1*a2^2*c9*c20*c22 + a1*a2^2*c10*c13*c28 + a1*a2^2*c10*c14*c27 + a1*a2^2*c10*c15*c26 + a1*a2^2*c10*c16*c25 - a1*a2^2*c10*c17*c24 - a1*a2^2*c10*c18*c23 - a1*a2^2*c10*c19*c22 - a1*a2^2*c10*c20*c21 + a1*a2^2*c11*c14*c26 - a1*a2^2*c11*c18*c22 + a1*a2^2*c12*c13*c26 + a1*a2^2*c12*c14*c25 - a1*a2^2*c12*c17*c22 - a1*a2^2*c12*c18*c21 + a1*a2*c1*c20 + a1*a2*c2*c19 - a1*a2*c2*c28 + a1*a2*c3*c18 + a1*a2*c4*c17 - a1*a2*c4*c26 - a1*a2*c5*c16 - a1*a2*c6*c15 + a1*a2*c6*c24 - a1*a2*c7*c14 - a1*a2*c8*c13 + a1*a2*c8*c22 + a1*a2*c9*c16*c28 - a1*a2*c9*c20*c24 + a1*a2*c10*c15*c28 + a1*a2*c10*c16*c27 - a1*a2*c10*c19*c24 - a1*a2*c10*c20*c23 + a1*a2*c11*c14*c28 + a1*a2*c11*c16*c26 - a1*a2*c11*c18*c24 - a1*a2*c11*c20*c22 + a1*a2*c12*c13*c28 + a1*a2*c12*c14*c27 + a1*a2*c12*c15*c26 + a1*a2*c12*c16*c25 - a1*a2*c12*c17*c24 - a1*a2*c12*c18*c23 - a1*a2*c12*c19*c22 - a1*a2*c12*c20*c21 + a1*c3*c20 + a1*c4*c19 - a1*c4*c28 - a1*c7*c16 - a1*c8*c15 + a1*c8*c24 + a1*c11*c16*c28 - a1*c11*c20*c24 + a1*c12*c15*c28 + a1*c12*c16*c27 - a1*c12*c19*c24 - a1*c12*c20*c23 + a2^3*c10*c14*c26 - a2^3*c10*c18*c22 + a2^2*c2*c18 - a2^2*c6*c14 + a2^2*c10*c14*c28 + a2^2*c10*c16*c26 - a2^2*c10*c18*c24 - a2^2*c10*c20*c22 + a2^2*c12*c14*c26 - a2^2*c12*c18*c22 + a2*c2*c20 + a2*c4*c18 - a2*c6*c16 - a2*c8*c14 + a2*c10*c16*c28 - a2*c10*c20*c24 + a2*c12*c14*c28 + a2*c12*c16*c26 - a2*c12*c18*c24 - a2*c12*c20*c22 + c4*c20 - c8*c16 + c12*c16*c28 - c12*c20*c24;
eqs(2) = a1*a2^2*c9 + a1*a2*c11 - a1*a2*c29 - a1*c31 + a2^2*c10 + a2*c12 - a2*c30 - c32;

end

%% The polynomial system
function cfg = retrieve_solver_cfg()
cfg = {};

% Number of coefficients, labelled as c1, c2, c3,...
cfg.numOfCoeff = 32;

% Number of variables, labeled as a1, a2, a3,...
cfg.numOfVars = 2;

% The index i of the selected variable, ai in the extra polynomial,
% ai - lambda.
% If set to -1, all variables will be tested one
% by one
cfg.hiddenVarNum = 2;

% (1) Either the size of polynomial combinations to be tested.
cfg.sizeOfCombs = [1;2];
% (2) Or the specific polynomial combination to be tested.
cfg.polyComb=[];
% (3) Or if both are given, the polycomb takes precedence over sizeofcombs.

% The number of rows to be GJ eliminated to obtained a reduced input
% polynomial system as an input to the generator.
cfg.noOfRowsToReduce = 0;

% The heuristic size of the template. There is no theoretical backing, to the best of our
% knowledge, governing the smallest template that can be generated.
% Hence one can start with a larger size and try to test by reducing the size
% of the template.
cfg.heurisiticTemplatesize = 117;
end