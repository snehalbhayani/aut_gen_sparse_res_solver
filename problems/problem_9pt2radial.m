function [ eqs, data0, eqs_data ] = problem_9pt2radial( data0 )

if nargin < 1 || isempty(data0)
    data0 = randi(50,72,1);
end

xx = create_vars(4);
f31 = xx(1);
f32 = xx(2);
l1 = xx(3);
l2 = xx(4);
f33 = 1;

vv = [ f31*l1 f32*l1 l1*l2 f31 f32 l1 l2 1]';

f11 = data0(1:8)'*vv;
f12 = data0(9:16)'*vv;
f13 = data0(17:24)'*vv;
f21 = data0(25:32)'*vv;
f22 = data0(33:40)'*vv;
f23 = data0(41:48)'*vv;
g5 = data0(49:56)'*vv;
g7 = data0(57:64)'*vv;
g9 = data0(65:72)'*vv;
F = [f11 f12 f13;f21 f22 f23;f31 f32 f33];
eqs = [l2*f13+g5;l2*f23+g7;f31*l1+g9;det(F)];
[C, M] = polynomials2matrix(eqs);
% eqs = rref(C) * M;

if nargout == 3
    xx = create_vars(4+72);
    data = xx(5:end);
    
f31 = xx(1);
f32 = xx(2);
l1 = xx(3);
l2 = xx(4);
f33 = 1;

vv = [ f31*l1 f32*l1 l1*l2 f31 f32 l1 l2 1]';

f11 = data(1:8)'*vv;
f12 = data(9:16)'*vv;
f13 = data(17:24)'*vv;
f21 = data(25:32)'*vv;
f22 = data(33:40)'*vv;
f23 = data(41:48)'*vv;
g5 = data(49:56)'*vv;
g7 = data(57:64)'*vv;
g9 = data(65:72)'*vv;
F = [f11 f12 f13;f21 f22 f23;f31 f32 f33];
eqs_data = [l2*f13+g5;l2*f23+g7;f31*l1+g9;det(F)];



end

