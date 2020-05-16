function[vars, hiddenvarnum, coeffconsts, sizeofcombs, polycomb, infinitePrec, eqs, theoreticalsolncnt, noofrowstoreduce, heuristictemplatesize] = problem_unsynch_relpose(data)
tic;
%%
numOfDataCoeff = 15*7;

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
for k = 1:7
    syms(strjoin({'a',num2str(k)},''));
    eval(strjoin({'xx(',num2str(k),') = ', 'a',num2str(k),';'},''));
end
data = transpose(data);
%%

%% Formatting the data structure

B = transpose([transpose(xx);data]);
p = mat2cell(B,1,ones(1,numel(B)));
if nargout >=7
    eqs = Eqs_problem_unsynch_relpose(p{:});
else
    eqs = [];
end
vars = transpose(xx);
if nargout >= 8
    theoreticalsolncnt = 16;
end


vars = transpose(xx);
coeffconsts = transpose(data);
hiddenvarnum = 7;
sizeofcombs=[2];
infinitePrec=2;
polycomb=[1;7];
noofrowstoreduce = 0;
heuristictemplatesize = 3600;
end