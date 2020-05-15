function[] = generate_five_pt_pose_solver(sizeOfBasisRed, folderName)
tic;
numOfDataCoeff = 36;

try
    rmdir(folderName, 's');
catch
    'No folder exists. Creating new one for now'
end

mkdir(folderName);
% 
% A simple wrapper for a five point pose estimation solver. 
% It simply reverse connects to a maple basis solver that 
% prints a matlab file 'FivePtCoeffMatricesName.m'.   

syms a1 a2 a3;

for k =1:numOfDataCoeff
    syms(strjoin({'c',num2str(k)},''));
    eval(strjoin({'data(',num2str(k),') = ', 'c',num2str(k),';'},''));
end

% data=randn(1,numOfDataCoeff);
% %==========================================================================
% %==========================================================================
% % This is the offline part where we generate the equations that are then
% % fed to maple.
% %==========================================================================
% %==========================================================================
xx = [a1;a2;a3];
E1 = reshape(data(1:9),3,3);
E2 = reshape(data(10:18),3,3);
E3 = reshape(data(19:27),3,3);
E4 = reshape(data(28:numOfDataCoeff),3,3);

E = a1 * E1 + a2 * E2+ a3 * E3 + E4;
% Creating and writing the essential matrix constraints to a txt file.
detEq = det(E);
eqs = [2*(E*transpose(E))*E - sum(diag(E*transpose(E)))*E];
eqs = [eqs(:);det(E)];
moreEqs='';
eqNameString = [];
for i = 1:length(eqs)
    constraintEq = char(expand(eqs(i)));
    eqNameString = [eqNameString ; string(strjoin({'eq',num2str(i)},''))];
    eqString = strjoin({' eq',num2str(i), ' := '},'');
    moreEqs = strjoin({moreEqs, eqString, constraintEq, ':'},'');
end
eqFile = fopen(strcat(folderName, '/Eqs.txt'),'w');
fprintf(eqFile ,'%s\n',moreEqs);
fprintf(eqFile, '%s\n', 'eqs := [',join(eqNameString,','),']:');
fprintf(eqFile ,'%s\n',strjoin({'vars := ',strrep(sym2str(xx),';',','),':'}));
fprintf(eqFile ,'%s\n',strjoin({'hiddenVarNumber := 3:'}));
fclose(eqFile);
disp("Printed the equations to a txt file ");
% Executing the maple code through a dos command. 
% We specify the folder name to be used, as that folder will house
% all of the data files that we will need.
% The needed data files include 1. the file containing matlab coefficient
% matrices for the polynomials to be solved. 2. a .txt file that houses all
% of the symbolic equations to be solved. 
baseMapleFile = fopen('MapleSparseBasisGeneratorTemplate.txt','r');
filetext = fileread('MapleSparseBasisGeneratorTemplate.txt');
fclose(baseMapleFile);

mapleFile = fopen('MapleSparseBasisGenerator.txt','w');
fprintf(mapleFile, '%s\n', strcat('solverFolderName := ',folderName, ':'));
fprintf(eqFile ,'%s\n',strjoin({'unHiddenVars := ',strrep(sym2str(xx(find(xx ~= a3))),';',','),' :'}));
fprintf(mapleFile, '%s\n', strcat('sizeOfBasisRed := ',num2str(sizeOfBasisRed),':'));
fprintf(mapleFile, '%s\n', strcat('sizeofcombs := [2]:'));
fprintf(mapleFile, '%s\n', strcat('numOfDataCoeff := ', num2str(numOfDataCoeff), ':'));
fprintf(mapleFile, '%s\n', strcat('varorder := [a1,a2,a3]:'));
fprintf(mapleFile, '%s\n', strcat('specificpolyorder:= []:'));

fprintf(mapleFile, '%s', filetext);
fclose(mapleFile);

if isunix
    mapleExecutor = strcat('maple  -q MapleSparseBasisGenerator.txt ');
elseif ispc
    mapleExecutor = strcat('cmaple  -q MapleSparseBasisGenerator.txt ');
else
    disp('Platform not supported')
end

dos(mapleExecutor);

%% 
% Using random data to estimate the rows and columns that can be removed to reduce the size of the monomial basis
mex reducemonbasis.c;
addpath(folderName);
data=randn(1,numOfDataCoeff);
% save('data','data');
% load('data');
% data = mat2cell(data, 1, ones(1,numOfDataCoeff));
[indicesToRemove, indxOfZeroCoeffMat, zeigvalindx, zeigvalindx2, solverTemplate, colpermutation,rowpermutation ] = reduce_mon_basis(data, numOfDataCoeff);

templatePart2 = fileread('PEP_solver_template.m');
templatePart1 = fileread(strcat(folderName, '/solver_template.m'));

solverFile = fopen(strcat(folderName, '/solver.m'),'w');
fprintf(solverFile, '%s \n', 'function[PEPsolutions] = solve(data, numofdatacoeff)');
fprintf(solverFile, '%s\n', templatePart1);
for i = 1:size(indicesToRemove,1)
    for j=1:size(indicesToRemove,2)
        fprintf(solverFile, '%s \n', strjoin({'indicesToRemove(', num2str(i), ',', num2str(j), ') = ', num2str(indicesToRemove(i,j)), ';'}, ''));
    end
end
fprintf(solverFile, '%s\n', strjoin({'indxOfZeroCoeffMat=',num2str(indxOfZeroCoeffMat), ';'},''));
fprintf(solverFile, '%s\n', strjoin({'colpermutation=',mat2str(colpermutation), ';'},''));
fprintf(solverFile, '%s\n', strjoin({'rowpermutation=',mat2str(rowpermutation), ';'},''));

fprintf(solverFile, '%s\n', strjoin({'zeigvalindx=',mat2str(zeigvalindx), ';'},''));
fprintf(solverFile, '%s\n', strjoin({'zeigvalindx2=',mat2str(zeigvalindx2), ';'},''));
fprintf(solverFile, '%s\n', strjoin({'solFromEigenVectors=transpose(',mat2str(solverTemplate), ');'},''));


fprintf(solverFile, '%s\n', templatePart2);
fclose(solverFile);

rmpath(folderName);
disp(toc);
end