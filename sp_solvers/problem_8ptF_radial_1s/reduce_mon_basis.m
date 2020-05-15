function[indicesToRemove, numOfRedCoeffMat, zeigvalindx, solverTemplate, extendedbasis, nullspacesize1, depdXcols1, nzrows1, depdCcols1, indepdCcols1, mdepdXcols1, mindepdXcols1, indepdCrows1, rowstorem1, nullspacesize2, depdXcols2, nzrows2, depdCcols2, indepdCcols2, mdepdXcols2, mindepdXcols2, indepdCrows2, rowstorem2, nullspacesize3, depdXcols3, nzrows3, depdCcols3, indepdCcols3, mdepdXcols3, mindepdXcols3, indepdCrows3, rowstorem3] = reducemonbasis(data, noofdatacoeff, foldername)
hiddenvarnumber = 3;
for i = 1:length(data)
eval(strjoin({'c',num2str(i),' = ', 'data(',num2str(i),');'},''));
end
t1 = c13*c25*c9-c17*c21*c9;
t2 = c11*c13*c25-c11*c17*c21+c13*c27*c9+c15*c25*c9-c17*c23*c9-c19*c21*c9-c1*c25+c21*c5;
t3 = c11*c13*c27+c11*c15*c25-c11*c17*c23-c11*c19*c21+c15*c27*c9-c19*c23*c9-c1*c27+c21*c7+c23*c5-c25*c3;
t4 = c11*c15*c27-c11*c19*c23+c23*c7-c27*c3;
t5 = c10*c13*c25-c10*c17*c21+c13*c26*c9+c14*c25*c9-c17*c22*c9-c18*c21*c9;
t6 = c10*c13*c27+c10*c15*c25-c10*c17*c23-c10*c19*c21+c11*c13*c26+c11*c14*c25-c11*c17*c22-c11*c18*c21+c12*c13*c25-c12*c17*c21+c13*c28*c9+c14*c27*c9+c15*c26*c9+c16*c25*c9-c17*c24*c9-c18*c23*c9-c19*c22*c9-c20*c21*c9+c1*c17-c1*c26-c13*c5-c2*c25+c21*c6+c22*c5;
t7 = c10*c15*c27-c10*c19*c23+c11*c13*c28+c11*c14*c27+c11*c15*c26+c11*c16*c25-c11*c17*c24-c11*c18*c23-c11*c19*c22-c11*c20*c21+c12*c13*c27+c12*c15*c25-c12*c17*c23-c12*c19*c21+c15*c28*c9+c16*c27*c9-c19*c24*c9-c20*c23*c9+c1*c19-c1*c28-c13*c7-c15*c5+c17*c3-c2*c27+c21*c8+c22*c7+c23*c6+c24*c5-c25*c4-c26*c3;
t8 = c11*c15*c28+c11*c16*c27-c11*c19*c24-c11*c20*c23+c12*c15*c27-c12*c19*c23-c15*c7+c19*c3+c23*c8+c24*c7-c27*c4-c28*c3;
t9 = c10*c13*c26+c10*c14*c25-c10*c17*c22-c10*c18*c21+c14*c26*c9-c18*c22*c9;
t10 = c10*c13*c28+c10*c14*c27+c10*c15*c26+c10*c16*c25-c10*c17*c24-c10*c18*c23-c10*c19*c22-c10*c20*c21+c11*c14*c26-c11*c18*c22+c12*c13*c26+c12*c14*c25-c12*c17*c22-c12*c18*c21+c14*c28*c9+c16*c26*c9-c18*c24*c9-c20*c22*c9+c1*c18-c13*c6-c14*c5+c17*c2-c2*c26+c22*c6;
t11 = c10*c15*c28+c10*c16*c27-c10*c19*c24-c10*c20*c23+c11*c14*c28+c11*c16*c26-c11*c18*c24-c11*c20*c22+c12*c13*c28+c12*c14*c27+c12*c15*c26+c12*c16*c25-c12*c17*c24-c12*c18*c23-c12*c19*c22-c12*c20*c21+c16*c28*c9-c20*c24*c9+c1*c20-c13*c8-c14*c7-c15*c6-c16*c5+c17*c4+c18*c3+c19*c2-c2*c28+c22*c8+c24*c6-c26*c4;
t12 = c11*c16*c28-c11*c20*c24+c12*c15*c28+c12*c16*c27-c12*c19*c24-c12*c20*c23-c15*c8-c16*c7+c19*c4+c20*c3+c24*c8-c28*c4;
t13 = c10*c14*c26-c10*c18*c22;
t14 = c10*c14*c28+c10*c16*c26-c10*c18*c24-c10*c20*c22+c12*c14*c26-c12*c18*c22-c14*c6+c18*c2;
t15 = c10*c16*c28-c10*c20*c24+c12*c14*c28+c12*c16*c26-c12*c18*c24-c12*c20*c22-c14*c8-c16*c6+c18*c4+c2*c20;
t16 = c12*c16*c28-c12*c20*c24-c16*c8+c20*c4;
t17 = c9;
t18 = c11-c29;
t19 = -c31;
t20 = c10;
t21 = c12-c30;
t22 = -c32;
M = zeros(3,17);
M(1,1) = t1;
M(1,2) = t2;
M(1,3) = t5;
M(1,4) = t3;
M(1,5) = t6;
M(1,6) = t9;
M(1,7) = t4;
M(1,8) = t7;
M(1,9) = t10;
M(1,10) = t13;
M(1,11) = t8;
M(1,12) = t11;
M(1,13) = t14;
M(1,14) = t12;
M(1,15) = t15;
M(1,17) = t16;
M(2,9) = t17;
M(2,12) = t18;
M(2,13) = t20;
M(2,14) = t19;
M(2,15) = t21;
M(2,17) = t22;
M(3,15) = 1;
M(3,16) = -1;
M = [rref(M(1:0,:)); M(1:end,:)];
ncinds = [1, 4, 7, 10, 13, 16, 19, 22, 25, 26, 28, 31, 34, 35, 37, 38, 40, 41, 43, 44, 49, 50];
for ncind = 1:22
    eval(strjoin({'nc', num2str(ncind), ' = ', 'M(ncinds(ncind));'}, '') )
end
ArrayOfCsNames = '';
sizeOfC = 1;
noOfVars = 3;
solFromEigenVectors = [s(1,1)(1) s(1,2)(1);];
Cs = zeros(1,1);
Cs(1,1) = Cred;
indicesToRemove = [];
indicesToSkip = fliplr(find(sum(transpose(solFromEigenVectors))~=0));
rowsToRemove = [];
colsToRemove = [];
rowsToRemove = rowsToRemove(1:max(find(rowsToRemove)), :);
colsToRemove = colsToRemove(1:max(find(colsToRemove)), :);
% if(size(rowsToRemove,1) > 0)
if 1 == 2
    indicesToRemove = [rowsToRemove,colsToRemove];
    
    noOfCoeffMatrices = size(Cs,2)/size(Cs,1);
    
    sizeOfReducedCs = sizeOfC - size(rowsToRemove,1);
    sizeOfM22 = size(rowsToRemove,1);
    numberOfCoeffMs = size(Cs,2)/size(Cs,1);
    maxReducedCoeffMs = 2 * numberOfCoeffMs - 1;
    colsToRemove = colsToRemove + ([1:numberOfCoeffMs] - 1) * sizeOfC;
    colsToRemove = colsToRemove(:);
    
    % Estimating M's
    M11 = Cs;
    M11(rowsToRemove,:) =[];
    M11(:,colsToRemove) =[];
    M22 = Cs(rowsToRemove, colsToRemove);
    M12 = Cs(:,colsToRemove);
    M12(rowsToRemove,:) = [];
    M21 = Cs(rowsToRemove,:);
    M21(:,colsToRemove) = [];
    reducedCs = [M11,zeros(sizeOfReducedCs,sizeOfReducedCs*(maxReducedCoeffMs-numberOfCoeffMs))];
    M22 = M22(1:sizeOfM22, 1:sizeOfM22);
    invM22 = inv(M22);
    indxOfZeroCoeffMat = [];
    for i=1:numberOfCoeffMs
        for j=1:numberOfCoeffMs
            coeffExp = i+j-1;
            temp = reducedCs(:,sizeOfReducedCs*(coeffExp-1) + 1:sizeOfReducedCs*coeffExp);
            tempM12 = M12(:,sizeOfM22*(i-1) + 1:sizeOfM22*i);
            tempM21 = M21(:,sizeOfReducedCs*(j-1) + 1:sizeOfReducedCs*j);
            
            buff = temp - (tempM12*invM22) * tempM21;
            reducedCs(:,sizeOfReducedCs*(coeffExp-1) + 1:sizeOfReducedCs*coeffExp) = buff;
            if( norm(buff,'fro') ~= 0)
                indxOfZeroCoeffMat(coeffExp) = 1;
            else
                indxOfZeroCoeffMat(coeffExp) = 0;
            end
        end
    end
    
    allReducedCs = [];
    for i = 1:maxReducedCoeffMs
        buff = reducedCs(:,sizeOfReducedCs*(i-1) + 1:sizeOfReducedCs*i);
        if(indxOfZeroCoeffMat(i) == 1)
            allReducedCs = [allReducedCs, buff];
            
        end
    end
    numOfRedCoeffMat = max(find(indxOfZeroCoeffMat));
    solForm = solFromEigenVectors;
    solForm(indicesToRemove(:,2),:) = [];
    sparseBasis(:,indicesToRemove(:,2)) = [];
    allCss = mat2cell(allReducedCs, sizeOfReducedCs, ones(1,numOfRedCoeffMat) * sizeOfReducedCs);
    
else
    sizeOfReducedCs = sizeOfC;
    numOfRedCoeffMat = size(Cs,2)/size(Cs,1);
    allCss = mat2cell(Cs, sizeOfReducedCs, ones(1,numOfRedCoeffMat ) * sizeOfReducedCs);
    indicesToRemove = [];
end


if rank(allCss{end}) < rank(allCss{1})
    allCss = fliplr(allCss);
end

A = []; B = [];
for i = 1:numOfRedCoeffMat -2
    tempA = []; tempB = [];
    for j = 1:numOfRedCoeffMat - 1
        if j == i + 1
            tempA = [tempA eye(sizeOfReducedCs)];
        else
            tempA = [tempA zeros(sizeOfReducedCs)];
        end
        if j == i
            tempB = [tempB eye(sizeOfReducedCs)];
        else
            tempB = [tempB zeros(sizeOfReducedCs)];
        end
    end
    A = [A;tempA];
    B = [B;tempB];
end

tempA = []; tempB = [];
for j = 1:numOfRedCoeffMat - 1
    tempA = [tempA -allCss{j}];
    if j == numOfRedCoeffMat -1
        tempB = [tempB allCss{j+1}];
    else
        tempB = [tempB zeros(sizeOfReducedCs)];
    end
end
A = [A; tempA];
B = [B; tempB];



extendedbasis = [sparseBasis; zeros(1,sizeOfReducedCs)];
for i = 1:numOfRedCoeffMat-2
    extendedbasis = [extendedbasis, [sparseBasis;ones(1,sizeOfReducedCs)*i]];
end


%% Removal of 0 eigen values via row algebra operations.

perms = [];
zeigvalindx=[];
colpermutation = [];
rowpermutation = [];

% remove_0_eigvalues;
% remove_more_0_eigvalues;


%% Generation of the template to extract solutions to individual
% unhidden variables from the eigen vector which will be a vector
% of monomials, which is an extension of the sparse bases.

% sizeoffinalres = size(extendedbasis,2);
% sizeoffinalres = min(sizeoffinalres, size(extendedbasis,2));
tempextendedbasis = extendedbasis(:,1:sizeoffinalres);
% tempextendedbasis = extendedbasis(:,end-sizeoffinalres+1:end);    
disp(tempextendedbasis);
numerr = 1e-10;

C0 = B; C1 = A;

A1 = -C0(1:end-sizeoffinalres, 1:end-sizeoffinalres);
A2 = -C0(1:end-sizeoffinalres, end-sizeoffinalres+1:end);
B1 = C1(end-sizeoffinalres+1:end, 1:end-sizeoffinalres);
B2 = C1(end-sizeoffinalres+1:end, end-sizeoffinalres+1:end);
X = B2 - (B1/A1) * A2;

[X, nullspacesize, depdXcols, nzrows, depdCcols, indepdCcols, mdepdXcols, mindepdXcols, indepdCrows, rowstorem] = remove_more_0_eigvalues_by_deflation(X);
nullspacesize1=nullspacesize;
depdXcols1=depdXcols;
nzrows1=nzrows;
depdCcols1=depdCcols;
indepdCcols1=indepdCcols;
mdepdXcols1=mdepdXcols;
mindepdXcols1=mindepdXcols;
indepdCrows1=indepdCrows;
rowstorem1=rowstorem;

[X, nullspacesize, depdXcols, nzrows, depdCcols, indepdCcols, mdepdXcols, mindepdXcols, indepdCrows, rowstorem] = remove_more_0_eigvalues_by_deflation(X);
nullspacesize2=nullspacesize;
depdXcols2=depdXcols;
nzrows2=nzrows;
depdCcols2=depdCcols;
indepdCcols2=indepdCcols;
mdepdXcols2=mdepdXcols;
mindepdXcols2=mindepdXcols;
indepdCrows2=indepdCrows;
rowstorem2=rowstorem;

[X, nullspacesize, depdXcols, nzrows, depdCcols, indepdCcols, mdepdXcols, mindepdXcols, indepdCrows, rowstorem] = remove_more_0_eigvalues_by_deflation(X);
nullspacesize3=nullspacesize;
depdXcols3=depdXcols;
nzrows3=nzrows;
depdCcols3=depdCcols;
indepdCcols3=indepdCcols;
mdepdXcols3=mdepdXcols;
mindepdXcols3=mindepdXcols;
indepdCrows3=indepdCrows;
rowstorem3=rowstorem;
    

pairs = transpose(combnk(setdiff(1:size(tempextendedbasis,2), rowstorem),2));

ops = eye(size(tempextendedbasis,1));
solverTemplate = zeros(size(tempextendedbasis));
for i = 1:size(tempextendedbasis,1) -1
    op = ops(:,i);
    for pair = pairs
        tmp = tempextendedbasis(:,pair)* [1;-1];
        if norm(abs(tmp) - op) == 0
            solverTemplate(:,pair(1)) = solverTemplate(:,pair(1)) + sum(tmp) * op;
            solverTemplate(:,pair(2)) = solverTemplate(:,pair(2)) -sum(tmp) * op;
            break;
        end
    end
end
%     Removoing the extended part of the sparseBasis
extendedbasis = tempextendedbasis;
%%


solverTemplate(size(extendedbasis,1),:) = [];

disp(" Parasitic eigenvalues have been removed. And size of final solver is ");
disp(size(extendedbasis));
end
