function [ eqs, data0, eqs_data ] = problem_relpose_7p_fr( data0 )

if nargin < 1 || isempty(data0)
 
    data0 = randi(50,56,1);
 
end

xx = create_vars(5);
f6 = xx(2);
f7 = xx(3);
f8 = xx(4);
l = xx(5);
w = xx(1);
f9 = 1;

h1 = data0(1:8);
h2 = data0(9:16);
h3 = data0(17:24);
h4 = data0(25:32);
h5 = data0(33:40);
h6 = data0(41:48);
h7 = data0(49:56);
vv = [l*f6 l*f7 l*f8  f6 f7 f8 l 1];
f1 = vv*h1;
f2 = vv*h2;
f3 = vv*h3;
f4 = vv*h4;
f5 = vv*h5;

F = [f1 f4 f7;f2 f5 f8;f3 f6 f9];
K = [1 0 0;0 1 0;0 0 w];
E = K'*F*K;

eqs = 2*(E*E')*E - sum(diag(E*E'))*E;
eqs = [eqs(:);det(F)];
eqs = [eqs;l*f3-vv*h6;l^2-vv*h7];

for iii = [1 2 4 5 9],
    cc = coeffs(eqs(iii));
    mm = monomials(eqs(iii));
    mm(1,:)=mm(1,:)/2;
    eqs(iii) = multipol(cc,mm);
end
for iii = [3 6 7 8],
    cc = coeffs(eqs(iii));
    mm = monomials(eqs(iii));
    mm(1,:) = mm(1,:)-1;
    mm(1,:)=mm(1,:)/2;
    eqs(iii) = multipol(cc,mm);
end



if nargout == 3
    xx = create_vars(5+56);
    data = xx(6:end);
    
    f6 = xx(2);
    f7 = xx(3);
    f8 = xx(4);
    l = xx(5);
    w = xx(1);
    f9 = 1;

h1 = data(1:8);
h2 = data(9:16);
h3 = data(17:24);
h4 = data(25:32);
h5 = data(33:40);
h6 = data(41:48);
h7 = data(49:56);
vv = [l*f6 l*f7 l*f8  f6 f7 f8 l 1];
f1 = vv*h1;
f2 = vv*h2;
f3 = vv*h3;
f4 = vv*h4;
f5 = vv*h5;

F = [f1 f4 f7;f2 f5 f8;f3 f6 f9];
K = [1 0 0;0 1 0;0 0 w];
E = K'*F*K;

eqs_data = 2*(E*E')*E - sum(diag(E*E'))*E;
eqs_data = [eqs_data(:);det(F)];
eqs_data = [eqs_data;l*f3-vv*h6;l^2-vv*h7];

for iii = [1 2 4 5 9],
    cc = coeffs(eqs_data(iii));
    mm = monomials(eqs_data(iii));
    mm(1,:)=mm(1,:)/2;
    eqs_data(iii) = multipol(cc,mm);
end
for iii = [3 6 7 8],
    cc = coeffs(eqs_data(iii));
    mm = monomials(eqs_data(iii));
    mm(1,:) = mm(1,:)-1;
    mm(1,:)=mm(1,:)/2;
    eqs_data(iii) = multipol(cc,mm);
end





end

