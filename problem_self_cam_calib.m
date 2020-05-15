function[eqs, vars, hiddenvarnum, coeffconsts, sizeofcombs, polycomb, infinitePrec] = problem_self_cam_calib(gen_instance)
tic;
%% Formatting the  structures -- coefficients and data 
numOfDataCoeff = 18;

if gen_instance == 1
    data0 = randn(1,numOfDataCoeff);
    data0(9) = 1;
    data0(18) = 1;
    
else
    for k = 1:numOfDataCoeff
        syms(strjoin({'c',num2str(k)},''));
        eval(strjoin({'data0(',num2str(k),') = ', 'c',num2str(k),';'},''));
    end
    eval(strjoin({'data0(9) = 1;'},''));
    eval(strjoin({'data0(18) = 1;'},''));
end

for k = 1:3
    syms(strjoin({'a',num2str(k)},''));
    eval(strjoin({'xx(',num2str(k),') = ', 'a',num2str(k),';'},''));
end
data0 = transpose(data0);

%% Formatting the data structure

nx = a1;
ny = a2;
l0 = a3;
l1 = a3;

% l0 = l1;
e = transpose([1 0 0]);
n = transpose([nx ny 1]); 
n = n/sqrt(nx^2 + ny^2 + 1);
a = cross(n,e);
b = cross(n,a);
K0 = diag([l0 l0 1]);
K1 = diag([l1 l1 1]);
K1inv = diag([1/l1 1/l1 1]);

for i = 1:9
    eval(strjoin({'H1(',num2str(i),') = ', 'data0(',num2str(i),');'},''));
end
for i = 10:18
    eval(strjoin({'H2(',num2str(i-9),') = ', 'data0(',num2str(i),');'},''));
end
H1 = reshape(H1, 3,3);
H2 = reshape(H2, 3,3);
at1 = K1inv  * H1 * K0 * a;
bt1 = K1inv  * H1 * K0 * b;
at2 = K1inv  * H2 * K1 * a;
bt2 = K1inv  * H2 * K1 * b;
eqs(1) = simplify(expand( (sqrt(nx^2 + ny^2 + 1))^3 * l1^2 * transpose(at1) * bt1));
eqs(2) = simplify(expand((nx^2 + ny^2 + 1)^2 * l1^2 * (transpose(at1) * at1 - transpose(bt1) * bt1)));
eqs(3) = simplify(expand( (sqrt(nx^2 + ny^2 + 1))^3 * l1^2 * transpose(at2) * bt2));

% [n,d] = numden(simplify(expand(transpose(at1) * bt1)));
% eqs(1) = n;
% [n,d] = numden(simplify(expand(transpose(at1) * at1 - transpose(bt1) * bt1)));
% eqs(2) = n;
% [n,d] = numden(simplify(expand(transpose(at2) * bt2)));
% eqs(3) = n;

vars = transpose(xx);

coeffconsts = transpose(data0);
hiddenvarnum = 2;
infinitePrec = 2;
sizeofcombs = [1;2];
polycomb=[];
end