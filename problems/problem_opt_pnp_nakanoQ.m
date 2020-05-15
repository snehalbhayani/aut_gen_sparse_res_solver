function [ eqs, data0, eqs_data ] = problem_opt_pnp_nakanoQ( data0 )

if nargin < 1 || isempty(data0)
    M = randi(200,9,9);
    M = M+M';
    data0 = M(:);
 
end

M = reshape(data0,9,9);
xx = create_vars(4);
a = xx(1);
b = xx(2);
c = xx(3);
d = xx(4);

R = [a^2+b^2-c^2-d^2 2*(b*c-a*d) 2*(b*d+a*c);...
    2*(b*c+a*d) a^2-b^2+c^2-d^2 2*(c*d-a*b);...
    2*(b*d-a*c) 2*(c*d+a*b) a^2-b^2-c^2+d^2];
r = R(:);
matMr = reshape(M*r,3,3);
P = R'*matMr-matMr'*R;
Q = matMr*R'-R*matMr';

eqs = [P(1,2);P(1,3);P(2,3);Q(1,2);Q(1,3);Q(2,3)];
eqs = [eqs;a^2+b^2+c^2+d^2-1];


if nargout == 3
    xx = create_vars(4+81);
    eqs_data = problem_opt_pnp_nakanoQ(xx(5:end));
end

