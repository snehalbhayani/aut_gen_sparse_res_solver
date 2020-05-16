function[vars, hiddenvarnum, coeffconsts, sizeofcombs, polycomb, infinitePrec, eqs, theoreticalsolncnt, noofrowstoreduce] = problem_opt_pnp_nakanoQ(data)
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

for k = 1:4
    syms(strjoin({'a',num2str(k)},''));
    eval(strjoin({'xx(',num2str(k),') = ', 'a',num2str(k),';'},''));
end
data = transpose(data);

%% Formatting the data structure

B = transpose([transpose(xx);data]);
p = mat2cell(B,1,ones(1,numel(B)));
if nargout >=7
    eqs = Eqs_opt_pnp_nakanoQ(p{:});
else
    eqs = [];
end

if nargout >= 8
    theoreticalsolncnt = 10;
end
vars = transpose(xx);

coeffconsts = transpose(data);
hiddenvarnum = -1;
infinitePrec = 2;
sizeofcombs = [3];
% polycomb=[2;6;7];
polycomb = [2;6;8];
noofrowstoreduce = size(eqs,2);
end