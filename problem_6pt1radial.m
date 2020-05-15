function[vars, hiddenvarnum, coeffconsts, sizeofcombs, polycomb, infinitePrec, eqs, theoreticalsolncnt, noofrowstoreduce, heurisitictemplatesize] = problem_6pt1radial(data)
tic;
%%
numOfDataCoeff = 54;

if nargin == 1
    if data == -1
        data = randn(1,numOfDataCoeff);
        for k = 1:4
            syms(strjoin({'a',num2str(k)},''));
            eval(strjoin({'xx(',num2str(k),') = ', 'a',num2str(k),';'},''));
        end
    else
        xx(1) = data(1);
        xx(2) = data(2);
        xx(3) = data(3);
        xx(4) = data(4);
        data = data(5:end);
    end
else
    for k = 1:numOfDataCoeff
        syms(strjoin({'c',num2str(k)},''));
        eval(strjoin({'data(',num2str(k),') = ', 'c',num2str(k),';'},''));
    end
        
    for k = 1:4
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
    eqs = Eqs_problem_6pt1radial(p{:});
else
    eqs = [];
end
vars = transpose(xx);
if nargout >= 8
    theoreticalsolncnt = 52;
end


vars = transpose(xx);
% vars = transpose([a3,a1,a2,a4]);
coeffconsts = transpose(data);
hiddenvarnum = 4;
sizeofcombs=[1];
infinitePrec=2;
polycomb = [2]; 
% polycomb = [5;12]; 
noofrowstoreduce = 0;
heurisitictemplatesize = 95;
end