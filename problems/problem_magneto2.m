function [ eqs, data0, eqs_data ] = problem_magneto2( data0 )


if nargin < 1 || isempty(data0)
   data0 = randi(50,3*4,1);
end

xx = create_vars(3);
x1 = data0(1:3);
x2 = data0(4:6);
y1 = data0(7:9);
y2 = data0(10:12);

Z = xx(1:3);
ZZ1 = Z-x1;
ZZ2 = Z-x2;
% Inverse matrices, i e m = MM1*y1 = MM2*y2
MM1 = 2*eye(3)*(-ZZ1'*ZZ1) + 3*(ZZ1*ZZ1');
MM2 = 2*eye(3)*(-ZZ2'*ZZ2) + 3*(ZZ2*ZZ2');
A = [MM1*y1 MM2*y2];
eqs(1,1) = det(A([2 3],:));
eqs(2,1) = det(A([1 3],:));
eqs(3,1) = det(A([1 2],:));
eqs =  [eqs;(MM1*y1)^2*(ZZ1'*ZZ1)-(MM2*y2)^2*(ZZ2'*ZZ2)];




if nargout == 3
    
    xx = create_vars(3+12);
    data = xx(4:end);
    
    x1 = data(1:3);
    x2 = data(4:6);
    y1 = data(7:9);
    y2 = data(10:12);

    Z = xx(1:3);
    ZZ1 = Z-x1;
    ZZ2 = Z-x2;
    % Inverse matrices, i e m = MM1*y1 = MM2*y2
    MM1 = 2*eye(3)*(-ZZ1'*ZZ1) + 3*(ZZ1*ZZ1');
    MM2 = 2*eye(3)*(-ZZ2'*ZZ2) + 3*(ZZ2*ZZ2');

    % Generate the 5 equations

    A = [MM1*y1 MM2*y2];
    eqs_data(1,1) = det(A([2 3],:));
    eqs_data(2,1) = det(A([1 3],:));
    eqs_data(3,1) = det(A([1 2],:));
    eqs_data =  [eqs_data;(MM1*y1)^2*(ZZ1'*ZZ1)-(MM2*y2)^2*(ZZ2'*ZZ2)];

    
    
end

