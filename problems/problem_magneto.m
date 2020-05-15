function [ eqs, data0, eqs_data ] = problem_magneto( data0 )


if nargin < 1 || isempty(data0)
   data0 = randi(50,3*4,1);
end

xx = create_vars(5);
x1 = data0(1:3);
x2 = data0(4:6);
y1 = data0(7:9);
y2 = data0(10:12);

Z = xx(1:3);
ZZ1 = Z-x1;
ZZ2 = Z-x2;
INZ1 = xx(4);
INZ2 = xx(5);
% Inverse matrices, i e m = MM1*y1 = MM2*y2
MM1 = 2*eye(3)*(-ZZ1'*ZZ1) + 3*(ZZ1*ZZ1');
MM2 = 2*eye(3)*(-ZZ2'*ZZ2) + 3*(ZZ2*ZZ2');

% Generate the 5 equations
eqs = INZ2*MM1*y1 - INZ1*MM2*y2;
eqs(4) = INZ1^2*(ZZ1'*ZZ1) - 1;
eqs(5) = INZ2^2*(ZZ2'*ZZ2) - 1;


if nargout == 3
    
    xx = create_vars(5+12);
    data = xx(6:end);
    
    x1 = data(1:3);
    x2 = data(4:6);
    y1 = data(7:9);
    y2 = data(10:12);

    Z = xx(1:3);
    ZZ1 = Z-x1;
    ZZ2 = Z-x2;
    INZ1 = xx(4);
    INZ2 = xx(5);
    % Inverse matrices, i e m = MM1*y1 = MM2*y2
    MM1 = 2*eye(3)*(-ZZ1'*ZZ1) + 3*(ZZ1*ZZ1');
    MM2 = 2*eye(3)*(-ZZ2'*ZZ2) + 3*(ZZ2*ZZ2');

    % Generate the 5 equations
    eqs_data = INZ2*MM1*y1 - INZ1*MM2*y2;
    eqs_data(4) = INZ1^2*(ZZ1'*ZZ1) - 1;
    eqs_data(5) = INZ2^2*(ZZ2'*ZZ2) - 1;


    
    
end

