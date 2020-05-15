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
