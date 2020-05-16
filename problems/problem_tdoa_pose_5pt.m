function[vars, hiddenvarnum, coeffconsts, sizeofcombs, polycomb, infinitePrec, eqs, actualsolno, noofrowstoreduce, degstotest] = problem_tdoa_pose_5pt(data)
tic;
%% Formatting the  structures -- coefficients and data 
numOfDataCoeff = 20;

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

for k = 1:5
    syms(strjoin({'a',num2str(k)},''));
    eval(strjoin({'xx(',num2str(k),') = ', 'a',num2str(k),';'},''));
end
data = transpose(data);
%% Formatting the data structure

B = transpose([transpose(xx);data]);
p = mat2cell(B,1,ones(1,numel(B)));
addpath('eqs');
if nargout >=7
    eqs = Eqs_tdoa_pose_5pt(p{:});
else
    eqs = [];
end
rmpath('eqs');

if nargout >= 8
    actualsolno = 20;
end
vars = transpose(xx);
coeffconsts = transpose(data);
hiddenvarnum = 4;
infinitePrec = 2;
sizeofcombs = [2];
polycomb=[];
noofrowstoreduce=size(eqs,2);
degstotest = [];
end