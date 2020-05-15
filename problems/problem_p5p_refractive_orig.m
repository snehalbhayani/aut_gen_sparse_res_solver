function [ eqs, data0, eqs_data ] = problem_p5p_refractive_orig( data0 )

if nargin < 1 || isempty(data0)
    data0 = randi(50,3*(5+5),1);
end
xx = create_vars(5);    
x = reshape(data0(1:3*5),3,5);
X = reshape(data0(3*5+1:end),3,5);


q = [1;xx(3:5)];
R = quat2rot(q);
C = [xx(1:2);0];

eqs = [];
for k = 1:5
    eqs = [eqs; cross(R'*x(:,k),X(:,k)-C)'*[0;0;1]];
end

if nargout == 3
    xx = create_vars(5+3*10);
    data = xx(6:end);
    x = reshape(data(1:3*5),3,5);
    X = reshape(data(3*5+1:end),3,5);

    q = [1;xx(3:5)];
    R = quat2rot(q);
    
    C = [xx(1:2);0];

    eqs_data = [];
    for k = 1:5
        eqs_data = [eqs_data; cross(R'*x(:,k),X(:,k)-C)'*[0;0;1]];
    end    
end

