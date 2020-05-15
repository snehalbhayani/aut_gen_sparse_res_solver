function[] = test_p2pr_refractive_solver()
% clc;
folder_name = 'sp_solvers/problem_p2pr_refractive';
problem_name = 'problem_p2pr_refractive';
addpath(folder_name);
addpath('eqs');
%% Tests are performed with the same value of R, T, F_true. All that changes is the sampled point correspondence measurements.
solverGenFunc = str2func(problem_name);

folderName = strcat('solvermats/',problem_name);
try
    rmdir(folderName, 's');
catch
    'No folder exists. Creating new one for now'
end
mkdir(folderName);

all_residual=[];
min_gt_error = [];
best_sols = [];
C0s = [];
C1s = [];
hiddenvarnumset = -1;
solverGenFunc = str2func(problem_name);
load("../scene_generator/problem_p2pr_refractive.mat");
load("../scene_generator/gt_p2pr.mat");

[varstemp, ~, coeffs, ~, ~, ~, symeqs, ~] = solverGenFunc();
iters = size(datas,2);
for index = 1:iters
    a1_t = g_t(1,index);
    a2_t = g_t(2,index);
    a3_t = g_t(3,index);
    a4_t = g_t(4,index);
    c = (1-a4_t ^2);
    s = 2* a4_t;
    R_t =  [[c,0,-s];[0,(1+a4_t^2),0];[s,0,c]];
    R_tnorm = norm(R_t);
    R_t = R_t/R_tnorm;
    
    data = datas(:,index);
    data = [data(1:3)/norm(data(1:3)); data(4:6)/norm(data(4:6)); data(7:16)];
    %     u = R_t'*reshape(data(1:6),3,2);
    %     X = reshape(data(7:12),3,2);
    %     X = transpose(R_t) * X;
    %     data = mat2cell([datas_p2pr(1:3,1);datas_p2pr(4:6,1); datas_p2pr(7:16,1)],ones(16,1),1)';
    %     disp(Eqs_problem_abspose_known_rot_refractive(a1_t, a2_t, a3_t, a4_t, data{:}));
    
    eqs = subs(symeqs, coeffs, transpose(data));
    if (mod(index,500)==0)
        fprintf("Tested %d iterations...",index);
    end
    residual = [];
    tic;
    [PEPSolutions, C0, C1, hiddenvarnum] = solver_problem_p2pr_refractive(data);
    C0s = [C0s, C0];
    C1s = [C1s, C1];
    elapsedTime = toc;
    if hiddenvarnumset == -1
        hiddenvarnumset = 1;
        hiddenvar = strjoin({'a',num2str(hiddenvarnum)},'');
        unhiddenvars = varstemp(find(varstemp ~= hiddenvar));
        vars = [unhiddenvars; hiddenvar];
    end

    t=[];r=[];
    PEPSolutions = transpose(PEPSolutions);
    for i = 1:size(PEPSolutions,2)
        sol = PEPSolutions([2,3,4,1],i);
        try
            t = [t, norm(abs([a1_t;a2_t;a3_t] -sol(1:3)))/norm(sol(1:3))];
            r = [r, norm(abs([a4_t] -sol(4)))/norm(sol(4))];
        catch
            disp("some error");
        end
    end
    sol = [r;t];
    try
        closestsol = PEPSolutions(:, find(min(sum(abs(sol))) == sum(abs(sol)),1));
%         temp = mat2cell([transpose(closestsol([2,3,4,1],:)),data'],1,1*ones(1,length(data)+length(closestsol)));
%         eqs = Eqs_problem_abspose_known_rot2(temp{:});
%         all_residual = [all_residual, max(abs(eqs))/norm(closestsol)];
        all_residual = [all_residual, max(abs(eval(subs(eqs, vars(1:4), closestsol([2,3,4,1])))))/norm(closestsol)];
        min_gt_error = [min_gt_error, [min(log10(sol(1,:)));min(log10(sol(2,:)))]];
        best_sols = [best_sols, closestsol];
        
        disp(mean(min_gt_error(1,:)));
        disp(mean(log10(all_residual)));
    catch
        disp("Some error");
    end
    
end
%%
save(strcat('sp_cvpr_', problem_name,'_C0s'), 'C0s', '-v7.3');
save(strcat('sp_cvpr_', problem_name,'_C1s'), 'C1s', '-v7.3');
save(strcat('sp_cvpr_',problem_name, '_best_sols'), 'best_sols');
save(strcat('sp_cvpr_',problem_name, '_min_gt_error'), 'min_gt_error');
save(strcat('sp_cvpr_',problem_name, '_all_residual'), 'all_residual');

rmpath(folder_name);
rmpath('eqs');
end