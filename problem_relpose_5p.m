function[vars, hiddenvarnum, coeffconsts, sizeofcombs, polycomb, infinitePrec, eqs, actualsolno, noofrowstoreduce] = problem_relpose_5p(data)
tic;
%% Formatting the data structures -- coefficients and data 
numOfDataCoeff = 36;

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
for k = 1:3
    syms(strjoin({'a',num2str(k)},''));
    eval(strjoin({'xx(',num2str(k),') = ', 'a',num2str(k),';'},''));
end
data = transpose(data);

%% Formatting the data structure
% Creating and writing the essential matrix constraints to a txt file.
if nargout >= 8
   actualsolno = 10;
end
E1 = reshape(data(1:9),3,3);
E2 = reshape(data(10:18),3,3);
E3 = reshape(data(19:27),3,3);
E4 = reshape(data(28:numOfDataCoeff),3,3);

E = a1 * E1 + a2 * E2+ a3 * E3 + E4;
% Creating and writing the essential matrix constraints to a txt file.
detEq = det(E);
eqs = [2*(E*transpose(E))*E - sum(diag(E*transpose(E)))*E];
eqs = [eqs(:);det(E)];

vars = transpose(xx);
coeffconsts = transpose(data);
hiddenvarnum = 1;
sizeofcombs = [2];
polycomb = [];
infinitePrec=2;
noofrowstoreduce = size(eqs,1);
end