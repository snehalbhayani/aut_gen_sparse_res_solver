function[] = generate_solver(problemName)
%% Step 1
folderName = "solvers/"+problemName;
try
    rmdir(folderName, 's');
catch
    'No folder exists. Creating new one for now'
end
mkdir(folderName);
solverGenFunc = str2func(problemName);
addpath("problems");
[eqsHandler, cfg] = solverGenFunc();
data = arrayfun(@(k) sym(char(strjoin({'c',num2str(k)},''))), [[1:cfg.numOfCoeff]], 'UniformOutput', false);
vars = arrayfun(@(k) sym(char(strjoin({'a',num2str(k)},''))), [[1:cfg.numOfVars]], 'UniformOutput', false);
p = [vars, data];
eqs = eqsHandler(p{:});
rmpath("problems");

%% Step 2
vars = [vars{:}];
hiddenvarnum = cfg.hiddenVarNum;
acthiddenvarnum = hiddenvarnum;
if hiddenvarnum ~= -1
    hiddenvar = strjoin({'a',num2str(hiddenvarnum)},'');
else
    hiddenvarnum = 1;
    hiddenvar = strjoin({'a',num2str(1)},'');
end
vars = transpose([strjoin({'a', num2str(hiddenvarnum)}, ''), vars(find(vars~=hiddenvar))]);
numOfDataCoeff = length(data);
eqsstring='';
eqnamestring = [];
for i = 1:length(eqs)
    constraintEq = char((eqs(i)));
    eqnamestring = [eqnamestring ; string(strjoin({'eq',num2str(i)},''))];
    eqString = strjoin({' eq',num2str(i), ' := '},'');
    eqsstring = strjoin({eqsstring, eqString, constraintEq, ':'},'');
end

%% Step 3
eqFile = fopen(strcat(folderName, '/Eqs.txt'),'w');
fprintf(eqFile ,'%s\n',eqsstring);
fprintf(eqFile, '%s\n', 'eqs := [',join(eqnamestring,','),']:');
try
    fprintf(eqFile ,'%s\n',strjoin({'vars := [',strjoin(arrayfun(@char, transpose(vars), 'uniform', 0),','), ']:' }, ''));
catch
    fprintf(eqFile ,'%s\n',strjoin({'vars := ',strrep(sym2str(vars), ';', ','), ':' }, ''));
end
fprintf(eqFile ,'%s\n',strjoin({'hiddenVarNumber := ',num2str(acthiddenvarnum),':'}));
fprintf(eqFile ,'%s\n',strjoin({'noofrowstoreduce := ',num2str(cfg.noOfRowsToReduce),':'}));
if acthiddenvarnum ~= -1
    
    try
        fprintf(eqFile ,'%s\n',strjoin({'unHiddenVars := [',strjoin(arrayfun(@char, transpose(vars(find(vars~=hiddenvar))), 'uniform', 0),','),']:'},''));
    catch
        fprintf(eqFile ,'%s\n',strjoin({'unHiddenVars := ',strrep(sym2str(vars(find(vars~=hiddenvar))), ';', ','), ':'},''));
    end
    
    
end
fclose(eqFile);
disp("Printed the equations to a txt file ");
%%
% Executing the maple code through a dos command.
% We specify the folder name to be used, as that folder will house
% all of the data files that we will need.
% The needed data files include 1. the file containing matlab coefficient
% matrices for the polynomials to be solved. 2. a .txt file that houses all
% of the symbolic equations to be solved.
baseMapleFile = fopen('MapleSparseBasisGeneratorTemplate2.txt','r');
filetext = fileread('MapleSparseBasisGeneratorTemplate2.txt');
fclose(baseMapleFile);

mapleFile = fopen(strcat(folderName, '/MapleSparseBasisGenerator.txt'),'w');
fprintf(mapleFile, '%s\n', strcat('solverFolderName := ',folderName, ':'));
if length(cfg.sizeOfCombs) == 1
    fprintf(mapleFile ,'%s\n',strjoin({'sizeofcombs := [',num2str(cfg.sizeOfCombs), ']:' }, ''));
else
    fprintf(mapleFile ,'%s\n',strjoin({'sizeofcombs := ',strrep(mat2str(cfg.sizeOfCombs),';',','), ':' }, ''));
end
fprintf(mapleFile ,'%s\n',strjoin({'heurisitictemplatesize := ',strrep(mat2str(cfg.heurisiticTemplatesize),';',','), ':' }, ''));
fprintf(mapleFile, '%s\n', strcat('noofdatacoeff := ',num2str(numOfDataCoeff), ':'));
fprintf(mapleFile, '%s\n', strcat('varorder := []:'));

if length(cfg.polyComb) == 1
    fprintf(mapleFile ,'%s\n',strjoin({'polycomb := [',num2str(cfg.polyComb), ']:' }, ''));
else
    fprintf(mapleFile ,'%s\n',strjoin({'polycomb := ',strrep(mat2str(cfg.polyComb),';',','), ':' }, ''));
end

fprintf(mapleFile, '%s', filetext);
fclose(mapleFile);

if isunix
    mapleExecutor = strjoin({'maple -q ', char(folderName), '/MapleSparseBasisGenerator.txt '}, '');
elseif ispc
    mapleExecutor = strjoin({'cmaple -q ', char(folderName), '/MapleSparseBasisGenerator.txt '}, '');
else
    disp('Platform not supported')
end

dos(mapleExecutor);
end