function[vars, hiddenvarnum, coeffconsts, sizeofcombs, polycomb, infinitePrec, eqs, theoreticalsolncnt, noofrowstoreduce, heurisitictemplatesize] = problem_9pt2radial(data)
tic;
%%
numOfDataCoeff = 72;

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
    eval(strjoin({'xxip(',num2str(k),') = ', 'a',num2str(k),';'},''));
end
data = transpose(data);
%%

%% Formatting the data structure

B = transpose([transpose(xxip);data]);
p = mat2cell(B,1,ones(1,numel(B)));
if nargout >=7
    eqs = Eqs_problem_9pt2radial(p{:});
else
    eqs = [];
end
vars = transpose(xxip);
if nargout >= 8
    theoreticalsolncnt = 24;
end

vars = transpose(xxip);
% vars = transpose([a3,a1,a2,a4]);
coeffconsts = transpose(data);
hiddenvarnum = 3;
sizeofcombs=[3];
infinitePrec=2;
polycomb=[1;2;3];
noofrowstoreduce = 0;
heurisitictemplatesize = 117;
end