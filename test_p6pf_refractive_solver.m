function[] = test_p6pf_refractive_solver()
% clc;
folder_name = 'sp_solvers/problem_p6pf_refractive';
problem_name = 'problem_p6pf_refractive';
addpath(folder_name);
addpath('eqs');
%% Tests are performed with the same value of R, T, F_true. All that changes is the sampled point correspondence measurements.

solverMatsFolderName = strcat('solvermats/',problem_name);
folderName = solverMatsFolderName;
try
    rmdir(folderName, 's');
catch
    'No folder exists. Creating new one for now'
end
mkdir(folderName);


% We record thr ground truth error for rot, translation and focal length.

all_residual=[];
min_gt_error = [];
best_sols = [];
C0s = [];
C1s = [];

hiddenvarnumset = -1;

load('../scene_generator/problem_p6pf_refractive.mat');
load('../scene_generator/gt_p6pf.mat');
solverGenFunc = str2func(problem_name);

[varstemp, hiddenvarnum, coeffs, sizeofcombs, polycomb, infinitePrec, symeqs, theoreticalsolncnt] = solverGenFunc();
iters = size(datas,2);
for index = 1:5000
    gt = g_t(:,index);
    t_gt = gt(5:7);
    R_gt = quat2rot(gt(1:4));
    t_gt = t_gt;
    gt(1:4) = gt(1:4)./gt(1);    
    data = transpose(datas(:,index));
    
    eqs = subs(symeqs, coeffs, data);
    tic;
    [PEPSolutions, C0, C1, hiddenvarnum] = solver(data);
    C0s = [C0s, C0];
    C1s = [C1s, C1];
    if hiddenvarnumset == -1
        hiddenvarnumset = 1;
        hiddenvar = strjoin({'a',num2str(hiddenvarnum)},'');
        unhiddenvars = varstemp(find(varstemp ~= hiddenvar));
        vars = [unhiddenvars; hiddenvar];
    end
    PEPSolutions = transpose(PEPSolutions);
    noofpepsolns = size(PEPSolutions,2);
    r=[];
    t=[];
    focal=[];
    for i = 1:noofpepsolns
        sol = PEPSolutions(:,i);
        try
            r = [r, norm(abs([1;abs(sol([1,2,3]))] - [abs(gt(1:4))]))/norm(abs(gt(1:4))) ];
            t = [t, norm(abs([sol([4,5])] - t_gt(1:2)))/norm(t_gt) ];
            focal = [focal, norm(abs(sol(6) - gt(9)))/norm(gt(9)) ];
        catch
            disp("Some error");
        end
    end
    sol = [r;t;focal];
    if (mod(index,1)==0)
        fprintf("Tested %d iterations...",index);
        try
            closestsol = PEPSolutions(:, find(min(sum(abs(sol))) == sum(abs(sol)),1));
            temp = mat2cell([transpose(closestsol),data],1,1*ones(1,length(data)+length(closestsol)));
            eqs = Eqs_red_problem_p6pf_refractive(temp{:});
            
            all_residual = [all_residual, max(abs(eqs))/norm(closestsol)];
            min_gt_error = [min_gt_error , [min(abs(sol(1,:))); min(abs(sol(2,:))); min(abs(sol(3,:)))]];
            best_sols = [best_sols, closestsol];
            disp(mean(log10(all_residual)));
            disp(mean(log10(min_gt_error(1,:) )));            
        catch
            disp("Some error");
        end
    end
    
end
save(strcat('sp_cvpr_', problem_name,'_C0s'), 'C0s', '-v7.3');
save(strcat('sp_cvpr_', problem_name,'_C1s'), 'C1s', '-v7.3');
save(strcat('sp_cvpr_',problem_name, '_best_sols'), 'best_sols');
save(strcat('sp_cvpr_',problem_name, '_min_gt_error'), 'min_gt_error');
save(strcat('sp_cvpr_',problem_name, '_all_residual'), 'all_residual');
rmpath(folder_name);
rmpath('eqs');
end