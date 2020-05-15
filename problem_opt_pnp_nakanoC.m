function[vars, hiddenvarnum, coeffconsts, sizeofcombs, polycomb, infinitePrec, eqs, theoreticalsolncnt, noofrowstoreduce, heurisitictemplatesize] = problem_opt_pnp_nakanoC(data)
tic;
%% Formatting the  structures -- coefficients and data 
numOfDataCoeff = 81;

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

for k = 1:3
    syms(strjoin({'a',num2str(k)},''));
    eval(strjoin({'xx(',num2str(k),') = ', 'a',num2str(k),';'},''));
end
data = transpose(data);

%% Formatting the data structure

B = transpose([transpose(xx);data]);
p = mat2cell(B,1,ones(1,numel(B)));

if nargout >=7
    eqs = Eqs_problem_opt_pnp_nakanoC(p{:});
else
    eqs = [];
end
vars = transpose(xx);
if nargout >= 8
    theoreticalsolncnt = 40;
end

coeffconsts = transpose(data);
hiddenvarnum = 3;
infinitePrec = 2;
sizeofcombs = [3];
polycomb=[4;6;8];
% polycomb=[1;2];
heurisitictemplatesize = 158;
noofrowstoreduce = 0;
% cs := [162, 17, 79, 152]: -- for opt_pnp_nakanoC
end