function[vars, hiddenvarnum, coeffconsts, sizeofcombs, polycomb, infinitePrec, eqs, actualsolno, noofrowstoreduce, heurisitictemplatesize ] = problem_relpose_7p_fr_1s_elim(data)
tic;
numOfDataCoeff = 60;

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
        data = data(length(xx)+1:end);
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

%% Formatting the data structure

B = transpose([transpose(xx);data]);
p = mat2cell(B,1,ones(1,numel(B)));
if nargout >=7
    eqs = Eqs_problem_relpose_7p_fr_1s_elim(p{:});
else
    eqs = [];
end
vars = transpose(xx);
if nargout >= 8
    actualsolno = 19;
end
coeffconsts = transpose(data);
hiddenvarnum = 4;
infinitePrec = 2;
sizeofcombs = [2];
polycomb=[1;2]; 
noofrowstoreduce = 0;
heurisitictemplatesize = 70;
end