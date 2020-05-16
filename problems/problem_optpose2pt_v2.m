function[vars, hiddenvarnum, coeffconsts, sizeofcombs, polycomb, infinitePrec, eqs, theoreticalsolncnt, noofrowstoreduce, d] = problem_optpose2pt_v2(data)
tic;
%%
numOfDataCoeff = 14;

if nargin == 1
    if data == -1
        data = randn(1,numOfDataCoeff);
        for k = 1:5
            syms(strjoin({'a',num2str(k)},''));
            eval(strjoin({'xx(',num2str(k),') = ', 'a',num2str(k),';'},''));
        end
    else
        xx(1) = data(1);
        xx(2) = data(2);
        xx(3) = data(3);
        xx(4) = data(4);
        xx(5) = data(5);
        data = data(6:end);
    end
else
    for k = 1:numOfDataCoeff
        syms(strjoin({'c',num2str(k)},''));
        eval(strjoin({'data(',num2str(k),') = ', 'c',num2str(k),';'},''));
    end
        
    for k = 1:5
        syms(strjoin({'a',num2str(k)},''));
        eval(strjoin({'xx(',num2str(k),') = ', 'a',num2str(k),';'},''));
    end
end
data = transpose(data);
%%

%% Formatting the data structure

B = transpose([transpose(xx);data]);
p = mat2cell(B,1,ones(1,numel(B)));
if nargout >=7
    eqs = Eqs_problem_optpose2pt_v2(p{:});
else
    eqs = [];
end
vars = transpose(xx);
if nargout >= 8
    theoreticalsolncnt = 24;
end

eqs = transpose(eqs);
vars = transpose(xx);
coeffconsts = transpose(data);
hiddenvarnum = 2;
sizeofcombs=[3];
infinitePrec=2;
polycomb=[1;2;3;7];
noofrowstoreduce =0;
% size(eqs,2);
d=200;
end