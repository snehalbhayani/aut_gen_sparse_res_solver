function [ eqs, data0, eqs_data ] = problem_p4p_fr( data0 )

if nargin < 1 || isempty(data0)
    data0 = randi(50,64,1);
end

xx = create_vars(4);
k = xx(1);
alfa = [xx(2:4);1];

M1 = reshape(data0(1:16),4,4);
M2 = reshape(data0(17:32),4,4);

M3 = reshape(data0(33:64),4,8);
P1 = M1*alfa;
P2 = M2*alfa;
v = [alfa(1:3);k*alfa(1:3);k;1];
P3 = M3*v;

eqs(1) = P1(1:3)'*P2(1:3);
eqs(2) = P1(1:3)'*P3(1:3);
eqs(3) = P2(1:3)'*P3(1:3);
eqs(4) = P1(1:3)'*P1(1:3)-P2(1:3)'*P2(1:3);

eqs = eqs(:);

if nargout == 3
    xx = create_vars(4+64);
    data = xx(5:end);

k = xx(1);
alfa = [xx(2:4);1];

M1 = reshape(data(1:16),4,4);
M2 = reshape(data(17:32),4,4);

M3 = reshape(data(33:64),4,8);
P1 = M1*alfa;
P2 = M2*alfa;
v = [alfa(1:3);k*alfa(1:3);k;1];
P3 = M3*v;

eqs_data(1) = P1(1:3)'*P2(1:3);
eqs_data(2) = P1(1:3)'*P3(1:3);
eqs_data(3) = P2(1:3)'*P3(1:3);
eqs_data(4) = P1(1:3)'*P1(1:3)-P2(1:3)'*P2(1:3);
eqs_data = eqs_data(:);


end

