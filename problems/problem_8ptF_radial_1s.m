function [ eqs, data0, eqs_data ] = problem_8ptF_radial_1s( data0 )

if nargin < 1 || isempty(data0)
    data0 = randi(50,4*8,1);
end

d1 = data0(1:4);
d2 = data0(4+(1:4));
d3 = data0(8+(1:4));
d4 = data0(12+(1:4));
d5 = data0(16+(1:4));
d6 = data0(20+(1:4));
d7 = data0(24+(1:4));
d8 = data0(28+(1:4));

xx = create_vars(2);

lam = xx(2);
f32 = xx(1);

vv = [lam*f32 lam f32 1]';
f11 = d1'*vv;
f21 = d2'*vv;
f31 = d3'*vv;
f12 = d4'*vv;
f22 = d5'*vv;
f13 = d6'*vv;
f23 = d7'*vv;
f33 = 1;
g2 = d8'*vv;

F = [f11 f12 f13;f21 f22 f23;f31 f32 f33];
eqs = det(F);
eqs = [eqs;lam*f31-g2];


if nargout == 3
    xx = create_vars(2+32);
    data = xx(3:end);
    

d1 = data(1:4);
d2 = data(4+(1:4));
d3 = data(8+(1:4));
d4 = data(12+(1:4));
d5 = data(16+(1:4));
d6 = data(20+(1:4));
d7 = data(24+(1:4));
d8 = data(28+(1:4));


lam = xx(2);
f32 = xx(1);

vv = [lam*f32 lam f32 1]';
f11 = d1'*vv;
f21 = d2'*vv;
f31 = d3'*vv;
f12 = d4'*vv;
f22 = d5'*vv;
f13 = d6'*vv;
f23 = d7'*vv;
f33 = 1;
g2 = d8'*vv;

F = [f11 f12 f13;f21 f22 f23;f31 f32 f33];
eqs_data = det(F);
eqs_data = [eqs_data;lam*f31-g2];


end

