function [ eqs, data0, eqs_data ] = problem_magneto3( data0 )


if nargin < 1 || isempty(data0)
   data0 = randi(50,3*4,1);
end

xx = create_vars(4);
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

eqs = MM1*y1-xx(4)*(MM2*y2);
eqs = [eqs;(ZZ2'*ZZ2)-xx(4)^2*(ZZ1'*ZZ1)];


if nargout == 3
    
    xx = create_vars(4+12);
    data = xx(5:end);
    
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

    eqs_data = MM1*y1-xx(4)*(MM2*y2);
    eqs_data = [eqs_data;(ZZ2'*ZZ2)-xx(4)^2*(ZZ1'*ZZ1)];

    
    
end

