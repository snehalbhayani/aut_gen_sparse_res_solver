function [eqs,data0,eqs_data] = problem_example(data0)

if nargin < 1 || isempty(data0)
    % no input, generate a random integer instance
    data0 = randi(10,3*4,1);
end

% Setup equation system
Q = reshape(data0,3,4);
xx = create_vars(3);
eqs = xx.^2 + Q * [xx;1];

% Setup equation with data as additional unknowns
if nargout == 3
    xx = create_vars(3+3*4);
    eqs_data = problem_example(xx(4:end));
end