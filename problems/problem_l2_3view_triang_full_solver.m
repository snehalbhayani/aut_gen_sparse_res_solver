function[eqs, vars, hiddenvarnum, coeffconsts, sizeofcombs, polycomb, infinitePrec] = problem_l2_3view_triang_full_solver(gen_instance)
tic;
%% Formatting the data structures -- coefficients and data 
numOfDataCoeff = 33;

if gen_instance == 1
    data = randn(1,numOfDataCoeff);
else
    for k = 1:numOfDataCoeff
        syms(strjoin({'c',num2str(k)},''));
        eval(strjoin({'data(',num2str(k),') = ', 'c',num2str(k),';'},''));
    end
end

for k = 1:9
    syms(strjoin({'a',num2str(k)},''));
    eval(strjoin({'xx(',num2str(k),') = ', 'a',num2str(k),';'},''));
end
data = transpose(data);

%% Formatting the data structure
for j = 10:42
    eval(strjoin({'x',num2str(j),' = ', 'c',num2str(j-9),';'},''));
    k = k+1;
end
for j = 1:9
    eval(strjoin({'x',num2str(j),' = ', 'a',num2str(j),';'},''));
    k = k+1;
end
eqs(1) =    2*x1 + x3*x7*x10 + x4*x7*x13 + x5*x9*x28 + x6*x9*x31 + x7*x16 + x9*x34 - 2*x37;
eqs(2) =    2*x2 + x3*x7*x11 + x4*x7*x14 + x5*x9*x29 + x6*x9*x32 + x7*x17 + x9*x35 - 2*x38;
eqs(3) =    x1*x7*x10 + x2*x7*x11 + 2*x3 + x5*x8*x19 + x6*x8*x22 + x7*x12 + x8*x25 - 2*x39;
eqs(4) =    x1*x7*x13 + x2*x7*x14 + 2*x4 + x5*x8*x20 + x6*x8*x23 + x7*x15 + x8*x26 - 2*x40;
eqs(5) =    x1*x9*x28 + x2*x9*x29 + x3*x8*x19 + x4*x8*x20 + 2*x5 + x8*x21 + x9*x30 - 2*x41;
eqs(6) =    x1*x9*x31 + x2*x9*x32 + x3*x8*x22 + x4*x8*x23 + 2*x6 + x8*x24 + x9*x33 - 2*x42;
eqs(7) =    x1*x3*x10 + x1*x4*x13 + x1*x16 + x2*x3*x11 + x2*x4*x14 + x2*x17 + x3*x12 + x4*x15 + x18;
eqs(8) =    x3*x5*x19 + x3*x6*x22 + x3*x25 + x4*x5*x20 + x4*x6*x23 + x4*x26 + x5*x21 + x6*x24 + x27;
eqs(9) =    x1*x5*x28 + x1*x6*x31 + x1*x34 + x2*x5*x29 + x2*x6*x32 + x2*x35 + x5*x30 + x6*x33 + x36;

vars = transpose(xx);
coeffconsts = transpose(data);
hiddenvarnum = 1;
infinitePrec = 2;
sizeofcombs = [4;5];
% polycomb = [1;2;6];
polycomb=[];
end