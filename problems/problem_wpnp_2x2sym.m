function [ eqs, data0, eqs_data ] = problem_wpnp_2x2sym( data0 )

if nargin < 1 || isempty(data0)
    data0 = randi(50,9,1);
end

x = create_vars(4);
VV = [0  -1i -1i 0;
     -1   0   0  1;
      1i  0   0  1i;
      0  -1   1  0];
%VV = VV/sqrt(2);
q = VV*x(1:4);

R = quat2rot(q);
A = diag(data0(1:3));
B = reshape(data0(4:end),2,3);
L = sum(sum((R(1:2,:)*A-B).^2));
eqs = diff(L);
eqs = eqs(:);

if nargout == 3
    xx = create_vars(4+9);
    data = xx(5:end);
    q = VV*x(1:4);
    R = quat2rot(q);
    A = diag(data(1:3));
    B = reshape(data(4:end),2,3);
    L = sum(sum((R(1:2,:)*A-B).^2));
    eqs_data = diff(L);
    eqs_data = eqs_data(1:4)';
end

