function[vars, hiddenvarnum, coeffconsts, sizeofcombs, polycomb, infinitePrec, eqs, actualsolno, noofrowstoreduce, heurisitictemplatesize ] = problem_relpose_7p_fr_1s_partial_elim(data)
tic;
%% Formatting the  structures -- coefficients and data 
numOfDataCoeff = 35;

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

B = transpose([transpose(xx);data]);
p = mat2cell(B,1,ones(1,numel(B)));
if nargout >=7
    eqs = Eqs_problem_relpose_7p_fr_1s_partial_elim(p{:});
else
    eqs = [];
end
vars = transpose(xx);
if nargout >= 8
    actualsolno = 19;
end
coeffconsts = transpose(data);
hiddenvarnum = 3;
infinitePrec = 2;
sizeofcombs = [2];
polycomb=[1;2]; 
noofrowstoreduce = 0;
heurisitictemplatesize = 41;
end