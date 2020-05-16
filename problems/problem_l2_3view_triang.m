function[vars, hiddenvarnum, coeffconsts, sizeofcombs, polycomb, infinitePrec, eqs, theoreticalsolncnt, noofrowstoreduce, degstotest] = problem_l2_3view_triang(data)
tic;
%% Formatting the data structures -- coefficients and data 
numOfDataCoeff = 24;

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

B = transpose([transpose(xx);data]);
p = mat2cell(B,1,ones(1,numel(B)));
if nargout >=7
    eqs = Eqs_problem_l2_3view_triang(p{:});
else
    eqs = [];
end
vars = transpose(xx);
if nargout >= 8
    theoreticalsolncnt = 40;
end

vars = transpose(xx);
coeffconsts = transpose(data);
hiddenvarnum = 7;
infinitePrec = 2;
sizeofcombs = [4];
polycomb = [];
% polycomb=[1;4;5;7;10];
noofrowstoreduce = 0;
degstotest = 10000;
end