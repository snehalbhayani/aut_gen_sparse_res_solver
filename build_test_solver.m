function[] = build_test_solver(build, test, itercnt)
%% Setup the solver generator
prompt = 'Enter the problem name...';
problemname = input(prompt); % absolute_pose_quivers_solver generate_optimal_PnPQ_solver 9pt2raddist_solver The name of the problem to be solved
foldername = strcat('sp_solvers/',problemname); % The name of the folder where will the generated solver files will be stored.
%% Build
fq = 1;
counter = 0;
while fq ~= 0 & counter < 500
    if build
        clearvars -except build fq counter foldername problemname itercnt test;
        generate_solver(strcat(foldername,''), problemname);
        fq = 0;
    end
    %% Test
    if test
        if nargin == 2
            itercnt = 10;
        end
        fprintf("Testing the solver..... with %d iterations.....\n", itercnt);
        fq = test_solver(problemname, foldername, itercnt);
    end
    counter = counter + 1;
end
end