function[vars, hiddenvarnum, coeffconsts, sizeofcombs, polycomb, infinitePrec, eqs, actualsolno, noofrowstoreduce, degstotest] = problem_optpose4pt_v2(data)
tic;
%% Formatting the  structures -- coefficients and data 
numOfDataCoeff = 28;

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
if nargout >= 7
    eqs = Eqs_optpose4pt_v2(p{:});
else
    eqs = [];
end
if nargout >= 8
   actualsolno = 28;
end

    vars = transpose(xx);

coeffconsts = transpose(data);
hiddenvarnum = 1;
infinitePrec = 2;
sizeofcombs = [3];
polycomb=[];
noofrowstoreduce = 0;
degstotest=[1;2;3];
end