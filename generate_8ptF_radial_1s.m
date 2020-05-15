function[] = generate_8ptF_radial_1s(sizeOfBasisRed, folderName)
tic;
try
    rmdir(folderName, 's');
catch
    'No folder exists. Creating new one for now'
end

mkdir(folderName);
numOfDataCoeff = 32;

syms a1 a2;


allConst = 'function [indicesToRemove] =  sw6ptPoseRadial(dummy';
for i = 1:numOfDataCoeff
    syms(strjoin({'c',num2str(i)},''));
    allConst = strjoin({allConst, ',', 'c',num2str(i)}, '');
    eval(strjoin({'data0(',num2str(i),') = ', 'c',num2str(i),';'},''));
end
allConst = strjoin({allConst,')'},'');
% data0 = transpose(data0);
% %==========================================================================
% %==========================================================================
% % This is the offline part where we generate the equations that are then
% % fed to maple.
% %==========================================================================
% %==========================================================================
% Creating and writing the essential matrix constraints to a txt file.

d1 = data0(1:4);
d2 = data0(4+(1:4));
d3 = data0(8+(1:4));
d4 = data0(12+(1:4));
d5 = data0(16+(1:4));
d6 = data0(20+(1:4));
d7 = data0(24+(1:4));
d8 = data0(28+(1:4));

xx = [a1;a2];

lam = xx(2);
f32 = xx(1);

vv = transpose([lam*f32 lam f32 1]);
f11 = d1*vv;
f21 = d2*vv;
f31 = d3*vv;
f12 = d4*vv;
f22 = d5*vv;
f13 = d6*vv;
f23 = d7*vv;
f33 = 1;
g2 = d8*vv;

F = [f11 f12 f13;f21 f22 f23;f31 f32 f33];
eqs = det(F);
eqs = [eqs;lam*f31-g2];


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
fprintf(eqFile ,'%s\n',strjoin({'hiddenVarNumber := 2:'}));
fclose(eqFile);disp("Printed the equations to a txt file ");
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
fprintf(mapleFile, '%s\n', strcat('paramsList := ',allConst,':'));
fprintf(eqFile ,'%s\n',strjoin({'unHiddenVars := [a1];'}));
fprintf(mapleFile, '%s\n', strcat('sizeOfBasisRed := ',num2str(sizeOfBasisRed),':'));
fprintf(mapleFile, '%s\n', strcat('sizeofcombs := [2,3]:'));
fprintf(mapleFile, '%s\n', strcat('numOfDataCoeff := ', num2str(numOfDataCoeff), ':'));

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
data = mat2cell(data, 1, ones(1,numOfDataCoeff));
% [indicesToRemove] = reduce_mon_basis(1, data{:});
    indxOfZeroCoeffMat = -1;
    indicesToRemove = [];
 
templatePart2 = fileread('PEP_solver_template.m');
templatePart1 = fileread(strcat(folderName, '/solver_template.m'));

solverFile = fopen(strcat(folderName, '/solver.m'),'w');
fprintf(solverFile, '%s', templatePart1);
for i = 1:size(indicesToRemove,1)
    for j=1:size(indicesToRemove,2)
        fprintf(solverFile, '%s \n', strjoin({'indicesToRemove(', num2str(i), ',', num2str(j), ') = ', num2str(indicesToRemove(i,j)), ';'}, ''));
    end
end
if indxOfZeroCoeffMat > -1
    fprintf(solverFile, '%s\n', strjoin({'indxOfZeroCoeffMat=',num2str(indxOfZeroCoeffMat), ';'},''));
end

fprintf(solverFile, '%s\n', templatePart2);
fclose(solverFile);

rmpath(folderName);%%
end