%% 
% Using random data to estimate the rows and columns that can be removed to reduce the size of the monomial basis
% mex reducemonbasis.c;
addpath(folderName);
data=randn(1,numOfDataCoeff);
[indicesToRemove, indxOfZeroCoeffMat, zeigvalindx, solverTemplate, extendedbasis] = reduce_mon_basis(data, numOfDataCoeff, folderName);

templatePart2 = fileread('PEP_solver_template.m');
templatePart1 = fileread(strcat(folderName, '/solver_template.m'));

%%
% A dynamically generated file for converting a given set of solutions to
% an eigenvector. This is then used to compare whether the given solution
% satisfies the constraint. 
subsmonfile = fopen(strcat(folderName, '/subsmon.m'),'w');
fprintf(subsmonfile, '%s\n', 'function [monsinstance] = subsmon(varvalues)');
fprintf(subsmonfile, '%s \n', strjoin({'extendedbasis=',mat2str(extendedbasis), ';'},''));
fprintf(subsmonfile, '%s \n', strjoin({'for r=1:size(extendedbasis,2)'},''));
fprintf(subsmonfile, '%s \n', strjoin({'temp=1;'},''));
fprintf(subsmonfile, '%s \n', strjoin({'for j=1:size(extendedbasis,1)'},''));
fprintf(subsmonfile, '%s \n', strjoin({'temp = temp * varvalues(j)^extendedbasis(j,r);'},''));
fprintf(subsmonfile, '%s \n', strjoin({'end'},''));
fprintf(subsmonfile, '%s \n', strjoin({'monsinstance(r,1) = temp;'},''));
fprintf(subsmonfile, '%s\n', 'end');
fprintf(subsmonfile, '%s\n', 'end');

%%
solverFile = fopen(strcat(folderName, '/solver.m'),'w');
fprintf(solverFile, '%s \n', 'function[PEPsolutions, C0, C1, hiddenvarnumber] = solve(data)');
fprintf(solverFile, '%s\n', templatePart1);
for i = 1:size(indicesToRemove,1)
    for j=1:size(indicesToRemove,2)
        fprintf(solverFile, '%s \n', strjoin({'indicesToRemove(', num2str(i), ',', num2str(j), ') = ', num2str(indicesToRemove(i,j)), ';'}, ''));
    end
end
classes = unique(extendedbasis(end,:));
monstocheck = extendedbasis(end,:);
fprintf(solverFile, '%s\n', strjoin({'classes =',mat2str(classes), ';'},''));
fprintf(solverFile, '%s \n', strjoin({'monstocheck=',mat2str(monstocheck), ';'},''));

fprintf(solverFile, '%s\n', strjoin({'indxOfZeroCoeffMat=',num2str(indxOfZeroCoeffMat), ';'},''));



fprintf(solverFile, '%s\n', strjoin({'zeigvalindx=',mat2str(zeigvalindx), ';'},''));
% fprintf(solverFile, '%s\n', strjoin({'premultiplier=',mat2str(premultiplier), ';'},''));
fprintf(solverFile, '%s\n', strjoin({'solFromEigenVectors=transpose(',mat2str(solverTemplate), ');'},''));

fprintf(solverFile, '%s\n', strjoin({'theoreticalsolncnt =',num2str(theoreticalsolncnt), ';'},''));
fprintf(solverFile, '%s\n', templatePart2);

fclose(solverFile);
%%
rmpath(folderName);
disp(toc);