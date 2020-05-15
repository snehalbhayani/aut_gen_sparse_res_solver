function [eqs,data0,eqs_data] = problem_example2_hv(data0)

zp = 30097;

if nargin < 1 || isempty(data0)
    % no input, generate a random integer instance
    data0 = randi(100,3*4*3,1);
end

% Setup equation system
Q1 = reshape(data0(1:12),3,4);
Q2 = reshape(data0(13:24),3,4);
Q3 = reshape(data0(25:36),3,4);

xx = create_vars(3);
eqs(1) = [xx(1:2);1]' * Q1 * [xx;1];
eqs(2) = [xx(1:2);1]' * Q2 * [xx;1];
eqs(3) = [xx(1:2);1]' * Q3 * [xx;1];

% Setup equation with data as additional unknowns
if nargout == 3
    xx0 = create_vars(3+3*4*3);
    xx = xx0(1:3);
    data = xx0(4:end);

    % Setup equation system
    Q1 = reshape(data(1:12),3,4);
    Q2 = reshape(data(13:24),3,4);
    Q3 = reshape(data(25:36),3,4);

    xx = create_vars(3);
    eqs_data(1) = [xx(1:2);1]' * Q1 * [xx;1];
    eqs_data(2) = [xx(1:2);1]' * Q2 * [xx;1];
    eqs_data(3) = [xx(1:2);1]' * Q3 * [xx;1];
end