function[vars, hiddenvarnum, coeffconsts, sizeofcombs, polycomb, infinitePrec, eqs, actualsolno, noofrowstoreduce,degstotest] = problem_homo_3pt(data)
tic;
numOfDataCoeff = 26;

if nargin == 1
    if data == -1
        data = randn(1,numOfDataCoeff);
    else
%         disp('Obtained data vector');
    end
else
    for k = 1:numOfDataCoeff
        syms(strjoin({'c',num2str(k)},''));
        eval(strjoin({'data(',num2str(k),') = ', 'c',num2str(k),';'},''));
    end
end


for k = 1:8
    syms(strjoin({'a',num2str(k)},''));
    eval(strjoin({'xx(',num2str(k),') = ', 'a',num2str(k),';'},''));
end
data = transpose(data);
%% Formatting the data structure

eqs(1) = a1*a6*a8*data(9)*data(24) - a1*a6*data(9)*data(21) + a2*a6*a8*data(12)*data(24) - a2*a6*data(12)*data(21) + a3*a6*a8*data(15)*data(24) - a3*a6*data(15)*data(21) - a4*data(12)*data(24) - a5*data(9)*data(24) + data(15)*data(21);
eqs(2) = -a1*a6*a7*data(9)*data(24) + a1*a6*data(9)*data(18) - a2*a6*a7*data(12)*data(24) + a2*a6*data(12)*data(18) - a3*a6*a7*data(15)*data(24) + a3*a6*data(15)*data(18) + a4*data(9)*data(24) - a5*data(12)*data(24) - data(15)*data(18);
eqs(3) = a1*a6*a8*data(10)*data(25) - a1*a6*data(10)*data(22) + a2*a6*a8*data(13)*data(25) - a2*a6*data(13)*data(22) + a3*a6*a8*data(16)*data(25) - a3*a6*data(16)*data(22) - a4*data(13)*data(25) - a5*data(10)*data(25) + data(16)*data(22);
eqs(4) = -a1*a6*a7*data(10)*data(25) + a1*a6*data(10)*data(19) - a2*a6*a7*data(13)*data(25) + a2*a6*data(13)*data(19) - a3*a6*a7*data(16)*data(25) + a3*a6*data(16)*data(19) + a4*data(10)*data(25) - a5*data(13)*data(25) - data(16)*data(19);
eqs(5) =  a1*a6*a8*data(11)*data(26) - a1*a6*data(11)*data(23) + a2*a6*a8*data(14)*data(26) - a2*a6*data(14)*data(23) + a3*a6*a8*data(17)*data(26) - a3*a6*data(17)*data(23) - a4*data(14)*data(26) - a5*data(11)*data(26) + data(17)*data(23);
eqs(6) = -a1*a6*a7*data(11)*data(26) + a1*a6*data(11)*data(20) - a2*a6*a7*data(14)*data(26) + a2*a6*data(14)*data(20) - a3*a6*a7*data(17)*data(26) + a3*a6*data(17)*data(20) + a4*data(11)*data(26) - a5*data(14)*data(26) - data(17)*data(20);
eqs(7) = a1^2 + a2^2 + a3^2 - 1;
eqs(8) = a4^2 + a5^2 - 1;

if nargout >= 8
    actualsolno = 8;
end

vars = transpose(xx);
coeffconsts = transpose(data);
hiddenvarnum = 8;
infinitePrec = 2;
sizeofcombs = [2;3];
polycomb=[];
degstotest = [];
noofrowstoreduce = 0;

end