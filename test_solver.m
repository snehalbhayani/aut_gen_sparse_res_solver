function fq = test_solver(problem_name, folder_name, iter_cnt)
% clc;
addpath(folder_name);
addpath('eqs');
%% Tests are performed with the same value of R, T, F_true. All that changes is the sampled point correspondence measurements.
threshold = 10^-3;
time_for_execution = [];
all_failling_residual=[];
all_results=[];
C0s = [];
C1s = [];
fq = 0;
hiddenvarnumset = -1;
load(strcat(problem_name));
datas = transpose(datas_p6p);
solverGenFunc = str2func(problem_name);
allsols = [];
[varstemp, hiddenvarnum, coeffs, sizeofcombs, polycomb, infinitePrec, symeqs, theoreticalsolncnt] = solverGenFunc();
hiddenvar = strjoin({'a',num2str(hiddenvarnum)},'');
varstemp = [strjoin({'a', num2str(hiddenvarnum)}, ''); varstemp(find(varstemp~=hiddenvar))];
vars = varstemp;
% autogenproblem = str2func(strjoin({char(problem_name), '_ag'}, '')); %problem_optpose4pt_v2_ag
% load('E:\sparsePEPSolver_v8\sp_solvers\problem_opt_pnp_nakanoC\sols.mat');
% datas = datas(:, randperm(5000));
for index = 1:iter_cnt
    data = transpose(datas(:,index));
%     data = randn(36,1)';
%         eqs = subs(symeqs, coeffs, data);
    
    res = [];
    [PEPSolutions, C0, C1, hiddenvarnum] = solver(data);
    C0s = [C0s ; C0];
    C1s = [C1s ; C1];
    PEPSolutions = transpose(PEPSolutions);
    for i = 1:size(PEPSolutions,2)
        sol = PEPSolutions(:,i); %[2,3,4,5,1]
        temp = mat2cell([transpose(sol),data],1,1*ones(1,length(data)+length(sol))); 
        eqs = Eqs_red_problem_p6pf_refractive(temp{:});        
        try
%             kthres = max(abs(eval(subs(eqs, vars, sol))))/norm(abs(sol));
            %             kthres = max(abs(evaluate(eqs, sol)))/norm(abs(sol));
                        kthres = max(abs(eqs))/norm(abs(sol));
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
    
    if isempty(res)
        all_failling_residual = [all_failling_residual, threshold+0.1];
    else
        all_failling_residual = [all_failling_residual, max(res)];
    end
    %             fprintf("%f   %f \n",mean(log10(all_results)), 100*length(find(log10(all_failling_residual) > 1e-3))/5000);
    
    
    if (mod(index-1,1)==0)
        fprintf("Tested %d iterations...",index);
        %                 if mean(log10(all_results)) > -9.4
        %                      fq = 1;
        %                      break;
        %                 end
        disp(mean(log10(all_results)));
    end
end
%%
disp(median(log10(all_results)));
disp(mean(log10(all_results)));

save(strcat(folder_name,'/C0s'), 'C0s', '-v7.3');
save(strcat(folder_name,'/C1s'), 'C1s', '-v7.3');
save(strcat(folder_name,'/residual'), 'all_results', '-v7.3');
save(strcat(folder_name,'/sols'), 'allsols', '-v7.3');

rmpath(folder_name);
fq;
end