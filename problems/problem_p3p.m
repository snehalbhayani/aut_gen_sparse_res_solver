function [ eqs, data0, eqs_data ] = problem_p3p( data0 )

if nargin < 1 || isempty(data0)
    data0 = randi(200,6*3,1);
end

xx = create_vars(3);
x = reshape(data0(1:9),3,3);
X = reshape(data0(10:18),3,3);
lambda = xx(1:3);
xl = x*diag(lambda);
d = sum((xl(:,[1 1 2])-xl(:,[2 3 3])).^2);
D = sum((X(:,[1 1 2])-X(:,[2 3 3])).^2);
eqs = (d-D)';

if nargout == 3
    xx = create_vars(3+6*3);
    eqs_data = problem_p3p(xx(4:end));
end

