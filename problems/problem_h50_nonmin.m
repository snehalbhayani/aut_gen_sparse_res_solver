function[eqsHandler, cfg] = problem_h50_nonmin()
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
function eqs = retrieve_eqs(a1,a2,a3,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30)
for k = 1:30
    eval(strjoin({'data(',num2str(k),') = ', 'c',num2str(k),';'},''));
end
H1 = reshape(data(1:10),10,1);
H2 = reshape(data(11:20),10,1);
H3 = reshape(data(21:30),10,1);

H = a3 * H1 + a2 * H2 + a1 * H3;

[h11, h21, h31, h12, h22, h32, h33, nx, ny, nz] = feval(@(x) x{:}, num2cell(H));
eqs = [h12 ^ 2 * nz ^ 2 + h22 ^ 2 * nz ^ 2 + h32 ^ 2 * nz ^ 2 - 2 * h32 * h33 * ny * nz + h33 ^ 2 * ny ^ 2 - ny ^ 2 - nz ^ 2 h11 * h12 * nz ^ 2 + h21 * h22 * nz ^ 2 + h31 * h32 * nz ^ 2 - h31 * h33 * ny * nz - h32 * h33 * nx * nz + h33 ^ 2 * nx * ny - nx * ny h11 * h12 * ny * nz - h12 ^ 2 * nx * nz + h21 * h22 * ny * nz - h22 ^ 2 * nx * nz + h31 * h32 * ny * nz - h31 * h33 * ny ^ 2 - h32 ^ 2 * nx * nz + h32 * h33 * nx * ny + nx * nz h11 ^ 2 * nz ^ 2 + h21 ^ 2 * nz ^ 2 + h31 ^ 2 * nz ^ 2 - 2 * h31 * h33 * nx * nz + h33 ^ 2 * nx ^ 2 - nx ^ 2 - nz ^ 2 h11 ^ 2 * ny * nz - h11 * h12 * nx * nz + h21 ^ 2 * ny * nz - h21 * h22 * nx * nz + h31 ^ 2 * ny * nz - h31 * h32 * nx * nz - h31 * h33 * nx * ny + h32 * h33 * nx ^ 2 - ny * nz h11 ^ 2 * ny ^ 2 - 2 * h11 * h12 * nx * ny + h12 ^ 2 * nx ^ 2 + h21 ^ 2 * ny ^ 2 - 2 * h21 * h22 * nx * ny + h22 ^ 2 * nx ^ 2 + h31 ^ 2 * ny ^ 2 - 2 * h31 * h32 * nx * ny + h32 ^ 2 * nx ^ 2 - nx ^ 2 - ny ^ 2 h11 * h12 * h32 * h33 * nz - h11 * h12 * h33 ^ 2 * ny - h12 ^ 2 * h31 * h33 * nz + h12 ^ 2 * h33 ^ 2 * nx + h21 * h22 * h32 * h33 * nz - h21 * h22 * h33 ^ 2 * ny - h22 ^ 2 * h31 * h33 * nz + h22 ^ 2 * h33 ^ 2 * nx + h11 * h12 * ny - h12 ^ 2 * nx + h21 * h22 * ny - h22 ^ 2 * nx + h31 * h32 * ny + h31 * h33 * nz - h32 ^ 2 * nx - h33 ^ 2 * nx + nx h11 ^ 2 * h32 * h33 * nz - h11 ^ 2 * h33 ^ 2 * ny - h11 * h12 * h31 * h33 * nz + h11 * h12 * h33 ^ 2 * nx + h21 ^ 2 * h32 * h33 * nz - h21 ^ 2 * h33 ^ 2 * ny - h21 * h22 * h31 * h33 * nz + h21 * h22 * h33 ^ 2 * nx + h11 ^ 2 * ny - h11 * h12 * nx + h21 ^ 2 * ny - h21 * h22 * nx + h31 ^ 2 * ny - h31 * h32 * nx - h32 * h33 * nz + h33 ^ 2 * ny - ny -h11 ^ 2 * h22 ^ 2 * nz - h11 ^ 2 * h32 ^ 2 * nz + h11 ^ 2 * h32 * h33 * ny + 2 * h11 * h12 * h21 * h22 * nz + 2 * h11 * h12 * h31 * h32 * nz - h11 * h12 * h31 * h33 * ny - h11 * h12 * h32 * h33 * nx - h12 ^ 2 * h21 ^ 2 * nz - h12 ^ 2 * h31 ^ 2 * nz + h12 ^ 2 * h31 * h33 * nx - h21 ^ 2 * h32 ^ 2 * nz + h21 ^ 2 * h32 * h33 * ny + 2 * h21 * h22 * h31 * h32 * nz - h21 * h22 * h31 * h33 * ny - h21 * h22 * h32 * h33 * nx - h22 ^ 2 * h31 ^ 2 * nz + h22 ^ 2 * h31 * h33 * nx + h11 ^ 2 * nz + h12 ^ 2 * nz + h21 ^ 2 * nz + h22 ^ 2 * nz + h31 ^ 2 * nz - h31 * h33 * nx + h32 ^ 2 * nz - h32 * h33 * ny - nz h11 ^ 2 * h22 ^ 2 * h33 ^ 2 - 2 * h11 * h12 * h21 * h22 * h33 ^ 2 + h12 ^ 2 * h21 ^ 2 * h33 ^ 2 - h11 ^ 2 * h22 ^ 2 - h11 ^ 2 * h32 ^ 2 - h11 ^ 2 * h33 ^ 2 + 2 * h11 * h12 * h21 * h22 + 2 * h11 * h12 * h31 * h32 - h12 ^ 2 * h21 ^ 2 - h12 ^ 2 * h31 ^ 2 - h12 ^ 2 * h33 ^ 2 - h21 ^ 2 * h32 ^ 2 - h21 ^ 2 * h33 ^ 2 + 2 * h21 * h22 * h31 * h32 - h22 ^ 2 * h31 ^ 2 - h22 ^ 2 * h33 ^ 2 + h11 ^ 2 + h12 ^ 2 + h21 ^ 2 + h22 ^ 2 + h31 ^ 2 + h32 ^ 2 + h33 ^ 2 - 1];

end

%% The polynomial system
function cfg = retrieve_solver_cfg()
cfg = {};
cfg.eqsInst = @retrieve_eqs;
% Number of coefficients, labelled as c1, c2, c3,...
cfg.numOfCoeff = 30;

% Number of variables, labeled as a1, a2, a3,...
cfg.numOfVars = 3;

% The index i of the selected variable, ai in the extra polynomial,
% ai - lambda.
% If set to -1, all variables will be tested one
% by one
cfg.hiddenVarNum = 3;

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
cfg.heurisiticTemplatesize = 100;
end