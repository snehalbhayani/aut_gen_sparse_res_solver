%% 
% Using random data to estimate the rows and columns that can be removed to reduce the size of the monomial basis
% mex reducemonbasis.c;
addpath(folderName);
data=randn(1,numOfDataCoeff);
[indicesToRemove, indxOfZeroCoeffMat, zeigvalindx, solverTemplate, extendedbasis, nullspacesize1, depdXcols1, nzrows1, depdCcols1, indepdCcols1, mdepdXcols1, mindepdXcols1, indepdCrows1, rowstorem1, nullspacesize2, depdXcols2, nzrows2, depdCcols2, indepdCcols2, mdepdXcols2, mindepdXcols2, indepdCrows2, rowstorem2, nullspacesize3, depdXcols3, nzrows3, depdCcols3, indepdCcols3, mdepdXcols3, mindepdXcols3, indepdCrows3, rowstorem3] = reduce_mon_basis(data, numOfDataCoeff, folderName);

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
fprintf(solverFile, '%s\n', strjoin({'nullspacesize1=',num2str(nullspacesize1), ';'},''));
fprintf(solverFile, '%s\n', strjoin({'nullspacesize2=',num2str(nullspacesize2), ';'},''));
fprintf(solverFile, '%s\n', strjoin({'nullspacesize3=',num2str(nullspacesize3), ';'},''));

fprintf(solverFile, '%s\n', strjoin({'depdXcols1=',mat2str(depdXcols1), ';'},''));
fprintf(solverFile, '%s\n', strjoin({'nzrows1=',mat2str(nzrows1), ';'},''));
fprintf(solverFile, '%s\n', strjoin({'depdCcols1=',mat2str(depdCcols1), ';'},''));
fprintf(solverFile, '%s\n', strjoin({'indepdCcols1=',mat2str(indepdCcols1), ';'},''));
fprintf(solverFile, '%s\n', strjoin({'mdepdXcols1=',mat2str(mdepdXcols1), ';'},''));
fprintf(solverFile, '%s\n', strjoin({'mindepdXcols1=',mat2str(mindepdXcols1), ';'},''));
fprintf(solverFile, '%s\n', strjoin({'indepdCrows1=',mat2str(indepdCrows1), ';'},''));
fprintf(solverFile, '%s\n', strjoin({'rowstorem1=',mat2str(rowstorem1), ';'},''));

fprintf(solverFile, '%s\n', strjoin({'depdXcols2=',mat2str(depdXcols2), ';'},''));
fprintf(solverFile, '%s\n', strjoin({'nzrows2=',mat2str(nzrows2), ';'},''));
fprintf(solverFile, '%s\n', strjoin({'depdCcols2=',mat2str(depdCcols2), ';'},''));
fprintf(solverFile, '%s\n', strjoin({'indepdCcols2=',mat2str(indepdCcols2), ';'},''));
fprintf(solverFile, '%s\n', strjoin({'mdepdXcols2=',mat2str(mdepdXcols2), ';'},''));
fprintf(solverFile, '%s\n', strjoin({'mindepdXcols2=',mat2str(mindepdXcols2), ';'},''));
fprintf(solverFile, '%s\n', strjoin({'indepdCrows2=',mat2str(indepdCrows2), ';'},''));
fprintf(solverFile, '%s\n', strjoin({'rowstorem2=',mat2str(rowstorem2), ';'},''));

fprintf(solverFile, '%s\n', strjoin({'depdXcols3=',mat2str(depdXcols3), ';'},''));
fprintf(solverFile, '%s\n', strjoin({'nzrows3=',mat2str(nzrows3), ';'},''));
fprintf(solverFile, '%s\n', strjoin({'depdCcols3=',mat2str(depdCcols3), ';'},''));
fprintf(solverFile, '%s\n', strjoin({'indepdCcols3=',mat2str(indepdCcols3), ';'},''));
fprintf(solverFile, '%s\n', strjoin({'mdepdXcols3=',mat2str(mdepdXcols3), ';'},''));
fprintf(solverFile, '%s\n', strjoin({'mindepdXcols3=',mat2str(mindepdXcols3), ';'},''));
fprintf(solverFile, '%s\n', strjoin({'indepdCrows3=',mat2str(indepdCrows3), ';'},''));
fprintf(solverFile, '%s\n', strjoin({'rowstorem3=',mat2str(rowstorem3), ';'},''));

fprintf(solverFile, '%s\n', strjoin({'zeigvalindx=',mat2str(zeigvalindx), ';'},''));
% fprintf(solverFile, '%s\n', strjoin({'premultiplier=',mat2str(premultiplier), ';'},''));
fprintf(solverFile, '%s\n', strjoin({'solFromEigenVectors=transpose(',mat2str(solverTemplate), ');'},''));

fprintf(solverFile, '%s\n', strjoin({'theoreticalsolncnt =',num2str(theoreticalsolncnt), ';'},''));
fprintf(solverFile, '%s\n', templatePart2);

fclose(solverFile);
%%
rmpath(folderName);
disp(toc);