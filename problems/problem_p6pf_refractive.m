function [ eqs, data0, eqs_data ] = problem_p6pf_refractive( data0 )

if nargin < 1 || isempty(data0)
    data0 = randi(50,3*(6+6),1);
end
    
x = reshape(data0(1:3*6),3,6);
X = reshape(data0(3*6+1:end),3,6);

xx = create_vars(6);
K = diag([1 1 xx(6)]);
q = [1;xx(1:3)];
R = quat2rot(q);
C = [xx(4:5);0];

eqs = [];
for k = 1:6
    eqs = [eqs; cross(R'*K*x(:,k),X(:,k)-C)'*[0;0;1]];
end



if nargout == 3
    xx = create_vars(6+3*12);
    data = xx(7:end);
    K = diag([1 1 xx(6)]);
    x = reshape(data(1:3*6),3,6);
    X = reshape(data(3*6+1:end),3,6);

    q = [1;xx(1:3)];
    R = quat2rot(q);
    
    C = [xx(4:5);0];

    eqs_data = [];
    for k = 1:6
        eqs_data = [eqs_data; cross(R'*K*x(:,k),X(:,k)-C)'*[0;0;1]];
    end
end

