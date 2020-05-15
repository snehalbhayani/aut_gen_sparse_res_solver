function[vars, hiddenvarnum, coeffconsts, sizeofcombs, polycomb, infinitePrec, eqs, actualsolno, noofrowstoreduce, degstotest] = problem_pose_35pt(data)
tic;
%% Formatting the data structures -- coefficients and data 
numOfDataCoeff = 19;

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
data = transpose(data);
for k = 1:2
    syms(strjoin({'a',num2str(k)},''));
    eval(strjoin({'xx(',num2str(k),') = ', 'a',num2str(k),';'},''));
end

%% Formatting the data structurey
B = transpose([transpose(xx);data]);
p = mat2cell(B,1,ones(1,numel(B)));
if nargout >= 7
    eqs = Eqs_pose_35pt(p{:});
else
    eqs = [];
end
if nargout >= 8
   actualsolno = 28;
end

    vars = transpose(xx);

coeffconsts = transpose(data);
hiddenvarnum = 2;
infinitePrec = 2;
sizeofcombs = [2];
polycomb=[];
noofrowstoreduce = size(eqs,2);
degstotest=[1;2;3];
end