function[vars, hiddenvarnum, coeffconsts, sizeofcombs, polycomb, infinitePrec, eqs, theoreticalsolncnt, noofrowstoreduce, heurisitictemplatesize] = problem_p6pf_refractive(data)
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
    for i = 1:numOfDataCoeff
        syms(strjoin({'c',num2str(i)},''));
        eval(strjoin({'data(',num2str(i),') = ', 'c',num2str(i),';'},''));
    end    
end

for k = 1:6
    syms(strjoin({'a',num2str(k)},''));
    eval(strjoin({'xx(',num2str(k),') = ', 'a',num2str(k),';'},''));
end
data = transpose(data);

%% Formatting the data structure
%% Formatting the data structure
% Creating and writing the essential matrix constraints to a txt file.
B = transpose([transpose(xx);data]);
p = mat2cell(B,1,ones(1,numel(B)));

if nargout >=7
    eqs = Eqs_red_problem_p6pf_refractive(p{:});
else
    eqs = [];
end
eqs = transpose(eqs);
vars = transpose(xx);

if nargout >= 8
    theoreticalsolncnt = 36;
end

vars = transpose(xx);
coeffconsts = transpose(data);
hiddenvarnum = 1;
sizeofcombs = [3];
polycomb=[1;2;3;4];
infinitePrec = 2;
noofrowstoreduce = 0;
heurisitictemplatesize = 300;
% size(eqs,2);
end