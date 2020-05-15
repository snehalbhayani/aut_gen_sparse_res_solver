function[] = generate_toy_solver(sizeOfBasisRed, folderName)
numOfDataCoeff = 35;
try
    rmdir(folderName, 's');
catch
    'No folder exists. Creating new one for now'
end
mkdir(folderName);
syms l1 l3 l2;

allConst = 'function [indicesToRemove] =  solve(dummy';
for i = 1:numOfDataCoeff
    syms(strjoin({'c',num2str(i)},''));
    allConst = strjoin({allConst, ',', 'c',num2str(i)}, '');
    eval(strjoin({'c(',num2str(i),') = ', 'c',num2str(i),';'},''));
end
allConst = strjoin({allConst,')'},'');


% eq1 = expand((l1^2-c1)*(l2^2-c2)*(l3^2-c3)*(l3^4-c4));
% eq2 = expand((l1^3-l1+9)*(l2^2-c5)*(-l1*l2*l3^2+l3^3+c10)*(l2*l3-l2-l3+1));
% eq3 = expand((-c6*l1^4+l1^2)*(-c7*l1^2+3*l1*l2+l2^2)*(l3^3-c8)*(l1*l2*l3+c9)*(l2*l3-l2-l3+1));
eq1 = expand((c21*l1-c1)*(c22*l2^2-c2)*(c23*l3-c3)*(c24*l3-c4)*(c25*l2-c5)*(c26*l1-c14)); 
eq2 = expand((c27*l1^2+c6)*(c28*l2-c7)*(c29*l2-c8)*(c30*l3^3-c9)*(-c17*l3+c31*l2)); 
eq3 = expand((c32*l1-c11)*(c33*l3-c12)*(c34*l2^2-c13)*(c35*l3-c20));
eqs = [eq1;eq2;eq3];

moreEqs='';
for i = 1:length(eqs)
    constraintEq = char(expand(eqs(i)));
    eqString = strjoin({' eq',num2str(i), ' := '},'');
    moreEqs = strjoin({moreEqs, eqString, constraintEq, ':'},'');
end
eqFile = fopen(strcat(folderName, '/Eqs.txt'),'w');
fprintf(eqFile ,'%s\n',moreEqs);
fprintf(eqFile, '%s\n', 'eqs := [eq1,eq2,eq3]:');
fprintf(eqFile ,'%s\n',strjoin({'vars := [l1,l2,l3] :'}));
fprintf(eqFile ,'%s\n',strjoin({'hiddenVar := l3:'}));
fclose(eqFile);



if isunix
    mapleExecutor = strcat('maple  -q auto_toy_sparse.txt ');
elseif ispc
    mapleExecutor = strcat('cmaple  -q auto_toy_sparse.txt ');
else
    disp('Platform not supported')
end

dos(mapleExecutor);

% %% 
% % Using random data to estimate the rows and columns that can be removed to reduce the size of the monomial basis
% mex reducemonbasis.c;
% addpath(folderName);
% data=randn(1,numOfDataCoeff);
% save('data','data');
% % load('data');
% data = mat2cell(data, 1, ones(1,numOfDataCoeff));
% [indicesToRemove] = reduce_mon_basis(1, data{:});
% indxOfZeroCoeffMat = indicesToRemove(1,3);
% indicesToRemove = indicesToRemove(:,1:2);
% 
% templatePart2 = fileread('PEP_solver_template.m');
% templatePart1 = fileread(strcat(folderName, '/solver_template.m'));
% 
% solverFile = fopen(strcat(folderName, '/solver.m'),'w');
% fprintf(solverFile, '%s', templatePart1);
% for i = 1:size(indicesToRemove,1)
%     for j=1:size(indicesToRemove,2)
%         fprintf(solverFile, '%s \n', strjoin({'indicesToRemove(', num2str(i), ',', num2str(j), ') = ', num2str(indicesToRemove(i,j)), ';'}, ''));
%     end
% end
% fprintf(solverFile, '%s\n', strjoin({'indxOfZeroCoeffMat=',num2str(indxOfZeroCoeffMat), ';'},''));
% 
% fprintf(solverFile, '%s\n', templatePart2);
% fclose(solverFile);
% 
% rmpath(folderName);
% disp(toc);
%%
end