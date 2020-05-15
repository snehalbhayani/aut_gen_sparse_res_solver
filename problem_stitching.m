function[vars, hiddenvarnum, coeffconsts, sizeofcombs, polycomb, infinitePrec, eqs, actualsolno, noofrowstoreduce, heuristictemplatesize] = problem_stitching2(data)
tic;
%% Formatting the  structures -- coefficients and data 
numOfDataCoeff = 32;
numofvars = 2;

if nargin == 1
    if data == -1
        data = randn(1,numOfDataCoeff);
        for k = 1:numofvars
            syms(strjoin({'a',num2str(k)},''));
            eval(strjoin({'xx(',num2str(k),') = ', 'a',num2str(k),';'},''));
        end
    else
        for k = 1:numofvars
            eval(strjoin({'xx(',num2str(k), ') = data(',num2str(k),');'}, ''));
        end
        data = data((numofvars+1):end);
    end
else
    for k = 1:numOfDataCoeff
        syms(strjoin({'c',num2str(k)},''));
        eval(strjoin({'data(',num2str(k),') = ', 'c',num2str(k),';'},''));
    end
    
    
    for k = 1:numofvars
        syms(strjoin({'a',num2str(k)},''));
        eval(strjoin({'xx(',num2str(k),') = ', 'a',num2str(k),';'},''));
    end
end
data = transpose(data);


%% Formatting the data structure

B = transpose([transpose(xx);data]);
p = mat2cell(B,1,ones(1,numel(B)));
addpath('eqs');
if nargout >=7
    eqs = Eqs_problem_stitching(p{:});
else
    eqs = [];
end
rmpath('eqs');

if nargout >= 8
    actualsolno = 18;
end
vars = transpose(xx);
coeffconsts = transpose(data);
hiddenvarnum = 1;
infinitePrec = 2;
sizeofcombs = [1];
polycomb=[1;2];
noofrowstoreduce=0;
heuristictemplatesize = 36;
end