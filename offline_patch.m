X = B1 - B2 * (A2 \ A1); 
[V,D] = eig(X);
EValues = diag(D); 
EVectors = V;
good = ~(isinf(EValues) | isnan(EValues));
EValues = EValues(good);
EVectors = EVectors(:,good);

% The eigen values and eigen vectors have been now extracted
PEPsolutions=[];
nonInfEValuesInd = ~isinf(EValues);
NinfEValues = EValues(nonInfEValuesInd);
NinfEVectors = EVectors(:,nonInfEValuesInd);
noOfEvalues = length(NinfEValues);
sizeOfEvectors = length(solForm);

% We basically, then iterate through all of the received eigenvalues and
% then try to remove those that do not satisfy the criterion for the
% corresponding eigenvectors to have a form that is the same as that of the
% monomial vector.
% In fact we also remove those eigenvalues and eigenvectors which give us
% solutions that have infinity value for atleast one variable.
g = 1;
% monstocheck(rowstorem1)=[];
% monstocheck(rowstorem2)=[];
allvarsextracted = sum(abs(solForm));
if length(find(allvarsextracted==0)) == 0
    
    for i = 1:noOfEvalues
        otherVarValues = NinfEVectors(:,i);
        for k = 1:noOfVars
            sols(i,k) = 1;
            for j = 1:sizeOfEvectors
                sols(i,k) = sols(i,k) * otherVarValues(j) ^ solForm(j,k);
            end
        end
        %     sols(i,1) = NinfEValues(i);
        sols(i,k + 1) = NinfEValues(i);
        PEPsolutions(g,:) = sols(i,:);
        g = g+1;
        
    end

else
    PEPsolutions = [];
end