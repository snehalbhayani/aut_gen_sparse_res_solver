function[eqsHandler, cfg] = problem_3view_calib_planar_scene()
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
% function eqs = retrieve_eqs(f1, nz, ny, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14, d15, d16, d17, d18, d19, d20, d21, d22, d23, d24, d25, d26, d27, d28, d29, d30, d31, d32, d33, d34, d35, d36, d37, d38, d39, d40, d41, d42, d43, d44, d45, d46, d47, d48, d49, d50, d51, d52, d53, d54, d55, d56, d57, d58, d59, d60, d61, d62, d63, d64, d65, d66, d67, d68, d69, d70, d71)
function eqs = retrieve_eqs(nz,ny,f1, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17, c18,c19)
% c4*c7 -- > c19
eqs = [-c4*c7*f1^3*ny^2*nz-c5*c8*f1^3*ny^2*nz-c7^2*f1^4*ny*nz-c8^2*f1^4*ny*nz+c1*c7*f1^3*ny+c2*c8*f1^3*ny+c4^2*f1^2*ny*nz+c4*c7*f1^3*nz+c5^2*f1^2*ny*nz+c5*c8*f1^3*nz-c6*c9*f1*ny^2*nz-c9^2*f1^2*ny*nz-c1*c4*f1^2-c2*c5*f1^2+c3*c9*f1*ny+c6^2*ny*nz+c6*c9*f1*nz-c3*c6, -c7^2*f1^4*ny^4*nz^2-c8^2*f1^4*ny^4*nz^2+c4^2*f1^2*ny^4*nz^2+4*c4*c7*f1^3*ny^3*nz^2+c5^2*f1^2*ny^4*nz^2+4*c5*c8*f1^3*ny^3*nz^2-c9^2*f1^2*ny^4*nz^2-2*c1*c4*f1^2*ny^3*nz-2*c1*c7*f1^3*ny^2*nz-2*c2*c5*f1^2*ny^3*nz-2*c2*c8*f1^3*ny^2*nz+4*c4*c7*f1^3*ny*nz^2+4*c5*c8*f1^3*ny*nz^2+c6^2*ny^4*nz^2+4*c6*c9*f1*ny^3*nz^2-c7^2*f1^4*ny^2+c7^2*f1^4*nz^2-c8^2*f1^4*ny^2+c8^2*f1^4*nz^2+c1^2*f1^2*ny^2-2*c1*c4*f1^2*ny*nz-2*c1*c7*f1^3*nz+c2^2*f1^2*ny^2-2*c2*c5*f1^2*ny*nz-2*c2*c8*f1^3*nz-2*c3*c6*ny^3*nz-2*c3*c9*f1*ny^2*nz-c4^2*f1^2*nz^2+2*c4*c7*f1^3*ny-c5^2*f1^2*nz^2+2*c5*c8*f1^3*ny+4*c6*c9*f1*ny*nz^2-c9^2*f1^2*ny^2+c9^2*f1^2*nz^2+c1^2*f1^2+c2^2*f1^2+c3^2*ny^2-2*c3*c6*ny*nz-2*c3*c9*f1*nz-c4^2*f1^2-c5^2*f1^2-c6^2*nz^2+2*c6*c9*f1*ny+c3^2-c6^2, -c13*c16*c18^2*f1^3*ny^6*nz^3-c14*c17*c18^2*f1^3*ny^6*nz^3+c15*c16^2*c18*f1^3*ny^6*nz^3+c15*c17^2*c18*f1^3*ny^6*nz^3+c10*c16*c18^2*f1^3*ny^5*nz^2+c11*c17*c18^2*f1^3*ny^5*nz^2-c12*c16^2*c18*f1^3*ny^5*nz^2-c12*c17^2*c18*f1^3*ny^5*nz^2-c13^2*c15*c18*f1*ny^6*nz^3+c13*c15^2*c16*f1*ny^6*nz^3-3*c13*c16*c18^2*f1^3*ny^4*nz^3-c14^2*c15*c18*f1*ny^6*nz^3+c14*c15^2*c17*f1*ny^6*nz^3-3*c14*c17*c18^2*f1^3*ny^4*nz^3+3*c15*c16^2*c18*f1^3*ny^4*nz^3+3*c15*c17^2*c18*f1^3*ny^4*nz^3+2*c10*c13*c15*c18*f1*ny^5*nz^2+c10*c13*c18^2*f1^2*ny^4*nz^2-c10*c15^2*c16*f1*ny^5*nz^2-2*c10*c15*c16*c18*f1^2*ny^4*nz^2+2*c10*c16*c18^2*f1^3*ny^3*nz^2+2*c11*c14*c15*c18*f1*ny^5*nz^2+c11*c14*c18^2*f1^2*ny^4*nz^2-c11*c15^2*c17*f1*ny^5*nz^2-2*c11*c15*c17*c18*f1^2*ny^4*nz^2+2*c11*c17*c18^2*f1^3*ny^3*nz^2+c12*c13^2*c18*f1*ny^5*nz^2-2*c12*c13*c15*c16*f1*ny^5*nz^2+2*c12*c13*c16*c18*f1^2*ny^4*nz^2+c12*c14^2*c18*f1*ny^5*nz^2-2*c12*c14*c15*c17*f1*ny^5*nz^2+2*c12*c14*c17*c18*f1^2*ny^4*nz^2-c12*c15*c16^2*f1^2*ny^4*nz^2-c12*c15*c17^2*f1^2*ny^4*nz^2-2*c12*c16^2*c18*f1^3*ny^3*nz^2-2*c12*c17^2*c18*f1^3*ny^3*nz^2-3*c13^2*c15*c18*f1*ny^4*nz^3+3*c13*c15^2*c16*f1*ny^4*nz^3-c13*c16*c18^2*f1^3*ny^4*nz-3*c13*c16*c18^2*f1^3*ny^2*nz^3-3*c14^2*c15*c18*f1*ny^4*nz^3+3*c14*c15^2*c17*f1*ny^4*nz^3-c14*c17*c18^2*f1^3*ny^4*nz-3*c14*c17*c18^2*f1^3*ny^2*nz^3+c15*c16^2*c18*f1^3*ny^4*nz+3*c15*c16^2*c18*f1^3*ny^2*nz^3+c15*c17^2*c18*f1^3*ny^4*nz+3*c15*c17^2*c18*f1^3*ny^2*nz^3-c10^2*c15*c18*f1*ny^4*nz-c10^2*c18^2*f1^2*ny^3*nz-2*c10*c12*c13*c18*f1*ny^4*nz+2*c10*c12*c15*c16*f1*ny^4*nz-c10*c13*c15^2*ny^4*nz^2+4*c10*c13*c15*c18*f1*ny^3*nz^2+2*c10*c13*c18^2*f1^2*ny^2*nz^2-2*c10*c15^2*c16*f1*ny^3*nz^2-4*c10*c15*c16*c18*f1^2*ny^2*nz^2+c10*c16*c18^2*f1^3*ny^3+c10*c16*c18^2*f1^3*ny*nz^2-c11^2*c15*c18*f1*ny^4*nz-c11^2*c18^2*f1^2*ny^3*nz-2*c11*c12*c14*c18*f1*ny^4*nz+2*c11*c12*c15*c17*f1*ny^4*nz-c11*c14*c15^2*ny^4*nz^2+4*c11*c14*c15*c18*f1*ny^3*nz^2+2*c11*c14*c18^2*f1^2*ny^2*nz^2-2*c11*c15^2*c17*f1*ny^3*nz^2-4*c11*c15*c17*c18*f1^2*ny^2*nz^2+c11*c17*c18^2*f1^3*ny^3+c11*c17*c18^2*f1^3*ny*nz^2+c12^2*c13*c16*f1*ny^4*nz+c12^2*c14*c17*f1*ny^4*nz+c12^2*c16^2*f1^2*ny^3*nz+c12^2*c17^2*f1^2*ny^3*nz+c12*c13^2*c15*ny^4*nz^2+2*c12*c13^2*c18*f1*ny^3*nz^2-4*c12*c13*c15*c16*f1*ny^3*nz^2+4*c12*c13*c16*c18*f1^2*ny^2*nz^2+c12*c14^2*c15*ny^4*nz^2+2*c12*c14^2*c18*f1*ny^3*nz^2-4*c12*c14*c15*c17*f1*ny^3*nz^2+4*c12*c14*c17*c18*f1^2*ny^2*nz^2-2*c12*c15*c16^2*f1^2*ny^2*nz^2-2*c12*c15*c17^2*f1^2*ny^2*nz^2-c12*c16^2*c18*f1^3*ny^3-c12*c16^2*c18*f1^3*ny*nz^2-c12*c17^2*c18*f1^3*ny^3-c12*c17^2*c18*f1^3*ny*nz^2-3*c13^2*c15*c18*f1*ny^2*nz^3+c13^2*c18^2*f1^2*ny^3*nz+3*c13*c15^2*c16*f1*ny^2*nz^3-c13*c16*c18^2*f1^3*ny^2*nz-c13*c16*c18^2*f1^3*nz^3-3*c14^2*c15*c18*f1*ny^2*nz^3+c14^2*c18^2*f1^2*ny^3*nz+3*c14*c15^2*c17*f1*ny^2*nz^3-c14*c17*c18^2*f1^3*ny^2*nz-c14*c17*c18^2*f1^3*nz^3-c15^2*c16^2*f1^2*ny^3*nz-c15^2*c17^2*f1^2*ny^3*nz+c15*c16^2*c18*f1^3*ny^2*nz+c15*c16^2*c18*f1^3*nz^3+c15*c17^2*c18*f1^3*ny^2*nz+c15*c17^2*c18*f1^3*nz^3+c10^2*c12*c18*f1*ny^3+c10^2*c15^2*ny^3*nz-c10^2*c18^2*f1^2*ny*nz-c10*c12^2*c16*f1*ny^3-4*c10*c12*c13*c18*f1*ny^2*nz+4*c10*c12*c15*c16*f1*ny^2*nz-2*c10*c13*c15^2*ny^2*nz^2+2*c10*c13*c15*c18*f1*ny*nz^2-c10*c13*c18^2*f1^2*ny^2+c10*c13*c18^2*f1^2*nz^2-c10*c15^2*c16*f1*ny*nz^2-2*c10*c15*c16*c18*f1^2*ny^2-2*c10*c15*c16*c18*f1^2*nz^2+c11^2*c12*c18*f1*ny^3+c11^2*c15^2*ny^3*nz-c11^2*c18^2*f1^2*ny*nz-c11*c12^2*c17*f1*ny^3-4*c11*c12*c14*c18*f1*ny^2*nz+4*c11*c12*c15*c17*f1*ny^2*nz-2*c11*c14*c15^2*ny^2*nz^2+2*c11*c14*c15*c18*f1*ny*nz^2-c11*c14*c18^2*f1^2*ny^2+c11*c14*c18^2*f1^2*nz^2-c11*c15^2*c17*f1*ny*nz^2-2*c11*c15*c17*c18*f1^2*ny^2-2*c11*c15*c17*c18*f1^2*nz^2-c12^2*c13^2*ny^3*nz-c12^2*c14^2*ny^3*nz+c12^2*c16^2*f1^2*ny*nz+c12^2*c17^2*f1^2*ny*nz+2*c12*c13^2*c15*ny^2*nz^2+c12*c13^2*c18*f1*ny*nz^2-2*c12*c13*c15*c16*f1*ny*nz^2+2*c12*c13*c16*c18*f1^2*ny^2+2*c12*c13*c16*c18*f1^2*nz^2+2*c12*c14^2*c15*ny^2*nz^2+c12*c14^2*c18*f1*ny*nz^2-2*c12*c14*c15*c17*f1*ny*nz^2+2*c12*c14*c17*c18*f1^2*ny^2+2*c12*c14*c17*c18*f1^2*nz^2+c12*c15*c16^2*f1^2*ny^2-c12*c15*c16^2*f1^2*nz^2+c12*c15*c17^2*f1^2*ny^2-c12*c15*c17^2*f1^2*nz^2-c13^2*c15*c18*f1*ny^2*nz-c13^2*c15*c18*f1*nz^3+c13^2*c18^2*f1^2*ny*nz+c13*c15^2*c16*f1*ny^2*nz+c13*c15^2*c16*f1*nz^3-c14^2*c15*c18*f1*ny^2*nz-c14^2*c15*c18*f1*nz^3+c14^2*c18^2*f1^2*ny*nz+c14*c15^2*c17*f1*ny^2*nz+c14*c15^2*c17*f1*nz^3-c15^2*c16^2*f1^2*ny*nz-c15^2*c17^2*f1^2*ny*nz-c10^2*c12*c15*ny^2+c10^2*c12*c18*f1*ny+c10^2*c15^2*ny*nz+c10^2*c15*c18*f1*nz+c10*c12^2*c13*ny^2-c10*c12^2*c16*f1*ny-2*c10*c12*c13*c18*f1*nz+2*c10*c12*c15*c16*f1*nz-c10*c13*c15^2*nz^2+2*c10*c13*c15*c18*f1*ny+c10*c15^2*c16*f1*ny-c11^2*c12*c15*ny^2+c11^2*c12*c18*f1*ny+c11^2*c15^2*ny*nz+c11^2*c15*c18*f1*nz+c11*c12^2*c14*ny^2-c11*c12^2*c17*f1*ny-2*c11*c12*c14*c18*f1*nz+2*c11*c12*c15*c17*f1*nz-c11*c14*c15^2*nz^2+2*c11*c14*c15*c18*f1*ny+c11*c15^2*c17*f1*ny-c12^2*c13^2*ny*nz-c12^2*c13*c16*f1*nz-c12^2*c14^2*ny*nz-c12^2*c14*c17*f1*nz+c12*c13^2*c15*nz^2-c12*c13^2*c18*f1*ny-2*c12*c13*c15*c16*f1*ny+c12*c14^2*c15*nz^2-c12*c14^2*c18*f1*ny-2*c12*c14*c15*c17*f1*ny-c13^2*c15*c18*f1*nz+c13*c15^2*c16*f1*nz-c14^2*c15*c18*f1*nz+c14*c15^2*c17*f1*nz-c10^2*c12*c15+c10*c12^2*c13-c10*c13*c15^2-c11^2*c12*c15+c11*c12^2*c14-c11*c14*c15^2+c12*c13^2*c15+c12*c14^2*c15];
% eqs = [d1*f1^4*ny*nz+d2*f1^3*ny^2*nz+d3*f1^3*ny+d4*f1^3*nz+d5*f1^2*ny*nz+d7*f1*ny^2*nz+d10*ny*nz+d6*f1^2+d8*f1*ny+d9*f1*nz+d11, f1^2*ny^2*d22+f1^3*nz*d19+f1^4*nz^2*d14+f1*ny*d29+f1^3*ny*d18+f1^4*ny^2*d13+nz*ny^3*d32+ny^4*nz^2*d31+f1*nz*d30+nz*ny*d34+f1^2*nz^2*d24+f1^4*ny^4*nz^2*d12+f1^3*ny^3*nz^2*d15+f1^3*ny^2*nz*d16+f1^3*ny*nz^2*d17+f1^2*ny^4*nz^2*d20+d36+f1^2*ny^3*nz*d21+ny*nz*f1^2*d23+f1*ny^3*nz^2*d26+f1*ny^2*nz*d27+f1*ny*nz^2*d28+f1^2*d25+ny^2*d33+nz^2*d35, nz*ny^3*d66+f1^2*nz^2*d52+f1*nz^3*d63+f1*ny*d62+f1^3*ny^3*d42+ny^2*nz^2*d67+f1*ny^3*d58+f1*nz*d64+f1^3*nz^3*d46+f1^2*ny^2*d50+nz*ny*d69+ny^4*nz^2*d65+f1*ny^4*nz^3*d55+f1*ny^4*nz*d56+f1*ny^3*nz^2*d57+d71+f1*ny^2*nz^3*d59+f1*ny^2*nz*d60+f1*ny*nz^2*d61+ny*nz*f1^2*d51+f1*ny^6*nz^3*d53+f1*ny^5*nz^2*d54+f1^3*ny^6*nz^3*d37+f1^3*ny^5*nz^2*d38+f1^3*ny^4*nz^3*d39+f1^3*ny^4*nz*d40+f1^3*ny^3*nz^2*d41+f1^3*ny^2*nz^3*d43+f1^3*ny^2*nz*d44+f1^3*ny*nz^2*d45+f1^2*ny^4*nz^2*d47+f1^2*ny^3*nz*d48+f1^2*ny^2*nz^2*d49+ny^2*d68+nz^2*d70];
end

%% The polynomial system
function cfg = retrieve_solver_cfg()
cfg = {};

% Number of coefficients, labelled as c1, c2, c3,...
cfg.numOfCoeff = 18;

% Number of variables, labeled as a1, a2, a3,...
cfg.numOfVars = 3;

% The index i of the selected variable, ai in the extra polynomial,
% ai - lambda.
% If set to -1, all variables will be tested one
% by one
cfg.hiddenVarNum = 2;

% (1) Either the size of polynomial combinations to be tested.
cfg.sizeOfCombs = [3];
% (2) Or the specific polynomial combination to be tested.
cfg.polyComb=[1;2;3];
% (3) Or if both are given, the polycomb takes precedence over sizeofcombs.

% The number of rows to be GJ eliminated to obtained a reduced input
% polynomial system as an input to the generator.
cfg.noOfRowsToReduce = 0;

% The heuristic size of the template. There is no theoretical backing, to the best of our
% knowledge, governing the smallest template that can be generated.
% Hence one can start with a larger size and try to test by reducing the size
% of the template.
cfg.heurisiticTemplatesize = 650;
end