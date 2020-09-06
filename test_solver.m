function [] = test_solver(problem_name, iter_cnt)
% clc;
addpath("solvers/"+problem_name);
addpath("problems");
%% Tests are performed with the same value of R, T, F_true. All that changes is the sampled point correspondence measurements.

all_results=[];
solverGenFunc = str2func(problem_name);
[eqsHandler, cfg] = solverGenFunc();

vars = arrayfun(@(k) sym(char(strjoin({'a',num2str(k)},''))), [[1:cfg.numOfVars]], 'UniformOutput', false);
data = arrayfun(@(k) sym(char(strjoin({'c',num2str(k)},''))), [[1:cfg.numOfCoeff]], 'UniformOutput', false);
p = [vars, data];
symeqs = eqsHandler(p{:});
hiddenvar = strjoin({'a',num2str(cfg.hiddenVarNum)},'');
vars = [vars{:}];
coeffs = [data{:}];
vars = transpose([strjoin({'a', num2str(cfg.hiddenVarNum)}, ''), vars(find(vars~=hiddenvar))]);
load('f.mat');
for index = 1:iter_cnt
%     data = randn(length(coeffs), 1);
    data = coeffdatas(index,:)';
%     data = data/norm(data);
    res = [];
    [PEPSolutions] = solver(data);
    eqs = subs(symeqs, coeffs, data');
    PEPSolutions = transpose(PEPSolutions);
    for i = 1:size(PEPSolutions,2)
        sol = PEPSolutions([1:end-1],i);
%         sol = [5;9/104;-16/24];
        try
            kthres = max(abs(eval(subs(eqs, vars, sol))))/norm(abs(sol));
            if kthres == 0
                kthres = 1e-20;
            end
            res = [res, kthres];
        catch
        end
    end
    all_results = [all_results, res];
    minres = min(log10(res));
    gt = PEPSolutions(:, log10(res) == minres);
    disp(min(log10(res)));
%     disp(gt);
end
%%
disp(mean(log10(all_results)));
rmpath("problems");
rmpath("solvers/"+problem_name);
end