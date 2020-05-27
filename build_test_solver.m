function[] = build_test_solver(build, test, itercnt)
%% Setup the solver generator
prompt = 'Enter the problem name...';
problemName = input(prompt, 's'); % absolute_pose_quivers_solver generate_optimal_PnPQ_solver 9pt2raddist_solver The name of the problem to be solved
%% Build
if build
    generate_solver(problemName);
end
%% Test
if test
    if nargin == 2
        itercnt = 10;
    end
    fprintf("Testing the solver..... with %d iterations.....\n", itercnt);
    test_solver(problemName, itercnt);
end
end