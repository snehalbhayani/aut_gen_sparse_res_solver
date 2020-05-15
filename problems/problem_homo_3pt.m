function [ eqs, data0, eqs_data ] = problem_homo_3pt( data0 )

if nargin < 1 || isempty(data0)
    data0 = randi(50,18,1);
end


x1 = data0(1:3);
y1 = data0(4:6);
w1 = data0(7:9);
x2 = data0(10:12);
y2 = data0(13:15);
w2 = data0(16:18);

xx = create_vars(8);
nx = xx(1);
ny = xx(2);
nz = xx(3);
cth = xx(4);
sth = xx(5);
d = xx(6);
tx = xx(7);
ty = xx(8);


eqs = [];
for iii = 1:3,
    eqs = [eqs;w1(iii)*y2(iii)-cth*w2(iii)*y1(iii)-sth*w2(iii)*x1(iii)-d*nz*w1(iii)*y2(iii)-d*nx*x1(iii)*y2(iii)-d*ny*y1(iii)*y2(iii)+d*nz*ty*w1(iii)*w2(iii)+d*nx*ty*w2(iii)*x1(iii)+d*ny*ty*w2(iii)*y1(iii)];
    eqs = [eqs;cth*w2(iii)*x1(iii)-w1(iii)*x2(iii)-sth*w2(iii)*y1(iii)+d*nz*w1(iii)*x2(iii)+d*nx*x1(iii)*x2(iii)+d*ny*x2(iii)*y1(iii)-d*nz*tx*w1(iii)*w2(iii)-d*nx*tx*w2(iii)*x1(iii)-d*ny*tx*w2(iii)*y1(iii)];
end
eqs = [eqs;nx^2+ny^2+nz^2-1];
eqs = [eqs;cth^2+sth^2-1];

if nargout == 3
    xx = create_vars(8+18);
    data = xx(9:end);
    

x1 = data(1:3);
y1 = data(4:6);
w1 = data(7:9);
x2 = data(10:12);
y2 = data(13:15);
w2 = data(16:18);

nx = xx(1);
ny = xx(2);
nz = xx(3);
cth = xx(4);
sth = xx(5);
d = xx(6);
tx = xx(7);
ty = xx(8);


eqs_data = [];
for iii = 1:3,
    eqs_data = [eqs_data;w1(iii)*y2(iii)-cth*w2(iii)*y1(iii)-sth*w2(iii)*x1(iii)-d*nz*w1(iii)*y2(iii)-d*nx*x1(iii)*y2(iii)-d*ny*y1(iii)*y2(iii)+d*nz*ty*w1(iii)*w2(iii)+d*nx*ty*w2(iii)*x1(iii)+d*ny*ty*w2(iii)*y1(iii)];
    eqs_data = [eqs_data;cth*w2(iii)*x1(iii)-w1(iii)*x2(iii)-sth*w2(iii)*y1(iii)+d*nz*w1(iii)*x2(iii)+d*nx*x1(iii)*x2(iii)+d*ny*x2(iii)*y1(iii)-d*nz*tx*w1(iii)*w2(iii)-d*nx*tx*w2(iii)*x1(iii)-d*ny*tx*w2(iii)*y1(iii)];
end
eqs_data = [eqs_data;nx^2+ny^2+nz^2-1];
eqs_data = [eqs_data;cth^2+sth^2-1];

end

