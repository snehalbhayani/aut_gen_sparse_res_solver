function[vars, hiddenvarnum, coeffconsts, sizeofcombs, polycomb, infinitePrec, eqs, theoreticalsolncnt, noofrowstoreduce, heurisitictemplatesize] = problem_pose_quiver(data)
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
    for k = 1:numOfDataCoeff
        syms(strjoin({'c',num2str(k)},''));
        eval(strjoin({'data(',num2str(k),') = ', 'c',num2str(k),';'},''));
    end
end

for k = 1:4
    syms(strjoin({'a',num2str(k)},''));
    eval(strjoin({'xx(',num2str(k),') = ', 'a',num2str(k),';'},''));
end
data = transpose(data);
%% Formatting the data structurey
B = transpose([transpose(xx);data]);
p = mat2cell(B,1,ones(1,numel(B)));

if nargout >=7
    eqs = Eqs_problem_pose_quiver(p{:});
else
    eqs = [];
end
eqs = transpose(eqs);
vars = transpose(xx);

if nargout >= 8
    theoreticalsolncnt = 20;
end
vars = transpose(xx);
coeffconsts =  transpose(data);
hiddenvarnum = 1;
infinitePrec = 2;
sizeofcombs = [3];
% polycomb = [2;3;4]; % sparse resultant size was 54 and eigen solver size was 34 with submatrix inversion by hiding variable a4.
polycomb=[1;2;4]; % with no submatrix inversion
noofrowstoreduce = 0;
heurisitictemplatesize = 850000;
end