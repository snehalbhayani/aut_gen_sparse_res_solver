function [] = test_solver(problem_name, iter_cnt)
clc;
addpath('solvers/'+problem_name);
addpath('problems');
%% Tests are performed with the same value of R, T, F_true. All that changes is the sampled point correspondence measurements.

all_failling_residual=[];
all_results=[];
load(strcat(problem_name));
solverGenFunc = str2func(problem_name);
allsols = [];
[varstemp, hiddenvarnum, coeffs, sizeofcombs, polycomb, infinitePrec, symeqs, theoreticalsolncnt] = solverGenFunc();
hiddenvar = strjoin({'a',num2str(hiddenvarnum)},'');
vars = [strjoin({'a', num2str(hiddenvarnum)}, ''); varstemp(find(varstemp~=hiddenvar))];
for index = 1:iter_cnt
    data = transpose(datas(:,index));
    
    res = [];
    [PEPSolutions] = solver(data);
    eqs = subs(symeqs, coeffs, data);
    PEPSolutions = transpose(PEPSolutions);
    for i = 1:size(PEPSolutions,2)
        sol = PEPSolutions([1:end-1],i);
        temp = mat2cell([transpose(sol),data],1,1*ones(1,length(data)+length(sol))); 
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

rmpath('problems');
rmpath('solvers/'+problem_name);
end