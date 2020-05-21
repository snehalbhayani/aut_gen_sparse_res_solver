function [] = test_solver(problem_name, iter_cnt)
% clc;
addpath('solvers/'+problem_name);
addpath('problems');
addpath('data');
%% Tests are performed with the same value of R, T, F_true. All that changes is the sampled point correspondence measurements.

all_failling_residual=[];
all_results=[];
load(strcat(problem_name));
solverGenFunc = str2func(problem_name);
allsols = [];
[eqsHandler, cfg] = solverGenFunc();

vars = arrayfun(@(k) sym(char(strjoin({'a',num2str(k)},''))), [[1:cfg.numOfVars]], 'UniformOutput', false);
data = arrayfun(@(k) sym(char(strjoin({'c',num2str(k)},''))), [[1:cfg.numOfCoeff]], 'UniformOutput', false);
p = [vars, data];
symeqs = eqsHandler(p{:});
hiddenvar = strjoin({'a',num2str(cfg.hiddenvarnum)},'');
vars = [vars{:}];
coeffs = [data{:}];
vars = transpose([strjoin({'a', num2str(cfg.hiddenvarnum)}, ''), vars(find(vars~=hiddenvar))]);

for index = 1:iter_cnt
    data = (datas(:,index));    
    res = [];
    [PEPSolutions] = solver(data);
    eqs = subs(symeqs, coeffs, data');
    PEPSolutions = transpose(PEPSolutions);
    for i = 1:size(PEPSolutions,2)
        sol = PEPSolutions([1:end-1],i);
        try
            kthres = max(abs(eval(subs(eqs, vars, sol))))/norm(abs(sol));
            if kthres == 0
                kthres = 1e-20;
            end
            res = [res, kthres];
        catch
            errors=1;
        end
    end
    allsols = [allsols, {PEPSolutions}];
    all_results = [all_results, res];
       
    disp(mean(log10(all_results)));
    
end
%%
disp(median(log10(all_results)));
disp(mean(log10(all_results)));
rmpath('data');
rmpath('problems');
rmpath('solvers/'+problem_name);
end