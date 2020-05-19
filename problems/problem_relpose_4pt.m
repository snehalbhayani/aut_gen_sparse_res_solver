function[vars, hiddenvarnum, coeffconsts, sizeofcombs, polycomb, infinitePrec, eqs, actualsolno, noofrowstoreduce, heurisitictemplatesize] = problem_relpose_4pt(data)
tic;
%% Formatting the  structures -- coefficients and data 
numOfDataCoeff = 18;

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
    eqs = Eqs_relpose_4pt(p{:});
    eqs = [eqs(1,1);eqs(1,2);eqs(1,3);eqs(1,4);eqs(1,5)];
else
    eqs = [];
end
rmpath('eqs');

if nargout >= 8
    actualsolno = 16;
end
vars = transpose(xx);
coeffconsts = transpose(data);
hiddenvarnum = 4;
infinitePrec = 2;
sizeofcombs = [1;2;3;4;5;6;7];
polycomb=[];
noofrowstoreduce=6;
heurisitictemplatesize = 249;
end