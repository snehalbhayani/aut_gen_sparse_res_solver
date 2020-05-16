function[vars, hiddenvarnum, coeffconsts, sizeofcombs, polycomb, infinitePrec, eqs, theoreticalsolncnt, noofrowstoreduce, d] = problem_accv_222(data)
tic;
%%
numOfDataCoeff = 36;

if nargin == 1
    if data == -1
        data = randn(1,numOfDataCoeff);
        for k = 1:3
            syms(strjoin({'a',num2str(k)},''));
            eval(strjoin({'xx(',num2str(k),') = ', 'a',num2str(k),';'},''));
        end
    else
        xx(1) = data(1);
        xx(2) = data(2);
        xx(3) = data(3);
        data = data(4:end);
    end
else
    for k = 1:numOfDataCoeff
        syms(strjoin({'c',num2str(k)},''));
        eval(strjoin({'data(',num2str(k),') = ', 'c',num2str(k),';'},''));
    end
        
    for k = 1:3
        syms(strjoin({'a',num2str(k)},''));
        eval(strjoin({'xx(',num2str(k),') = ', 'a',num2str(k),';'},''));
    end
end
data = transpose(data);
%%

B = transpose([transpose(xx);data]);
p = mat2cell(B,1,ones(1,numel(B)));
if nargout >=7
    eqs = Eqs_problem_accv_222(p{:});
else
    eqs = [];
end
vars = transpose(xx);
if nargout >= 8
    theoreticalsolncnt = 54;
end

eqs = transpose(eqs);
vars = transpose(xx);
coeffconsts = transpose(data);
hiddenvarnum = 2;
sizeofcombs=[2];
infinitePrec=2;
polycomb=[];
noofrowstoreduce =0;
% size(eqs,2);
d=1000000;
end