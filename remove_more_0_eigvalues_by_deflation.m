function[X, nullspacesize, depdXcols, nzrows, depdCcols, indepdCcols, mdepdXcols, mindepdXcols, indepdCrows, rowstorem] = remove_more_0_eigvalues_by_deflation(X)
%%
% returning parameter
% nullspacesize
% depdXcols
% nzrows
% depdCcols
% indepdCcols
% mdepdXcols
% mindepdXcols
% depdCrows
% rowstorem
numerr = 1e-10;
% nullspacesize = size(X,2) - rank(X);
nullspacesize = 0;
if nullspacesize ~= 0
    while true
        ps=zeros(1,size(X,2));
        for p = randperm(size(X,2))
            t = X; t(:,p)= [];
            ps(p) = norm(t*pinv(t)*X(:,p) - X(:,p));
        end
        depdXcols = find(ps/norm(ps) < 1e-10);
        %     depdXcols = find(ps < 1e-10);
        for p = setdiff(randperm(size(X,2)), depdXcols)
            if rank(X(:,[depdXcols,p])) == rank(X) && rank(X(:,[depdXcols,p])) == rank(X(:,depdXcols))
                break;
            end
        end
        
        depdXcols = [depdXcols,p];
        C = X(:,depdXcols);
        nzrows = find(sum(abs(C)') > 1e-10);
%         nzrows = 1:size(C,1);
        C = C(nzrows,:);
        
        ranktocheck = rank(C);
        indepdCcols = randperm(size(C,2),ranktocheck);
        while(rank(C(:, indepdCcols)) ~= ranktocheck)
            indepdCcols = randperm(size(C,2),ranktocheck);
        end
        
        depdCcols = setdiff(1:size(C,2), indepdCcols);
        mdepdXcols = depdXcols(depdCcols);
        mindepdXcols = depdXcols(indepdCcols);
        
        Csliced = C(:,indepdCcols);
        ranktocheck = rank(Csliced);
        indepdCrows = randperm(size(Csliced,1),ranktocheck);
        for dummy = 1:500
            indepdCrows = randperm(size(Csliced,1),ranktocheck);
            v = null(C(indepdCrows,:));
            Xp1 = zeros(size(X,2), nullspacesize);
            Xp1(mdepdXcols,:) = -eye(nullspacesize);
            Xp1(mindepdXcols,:) = (C(indepdCrows,indepdCcols)\C(indepdCrows,depdCcols));
            disp('---------------------------------------');
            if norm(X*Xp1,'fro') < numerr
                break;
            end
            disp('--------------------------------------- end');
        end
        if norm(X*Xp1,'fro') < numerr
            break;
        end
    end
    rowstorem = mdepdXcols;
    Xp2 = X(rowstorem,:);
    X = X + Xp1 * Xp2;
    X(mdepdXcols,:)=[];
    X(:,mdepdXcols)=[];
else
    nullspacesize=0;
    depdXcols=[];
    nzrows=[];
    depdCcols=[];
    indepdCcols=[];
    mdepdXcols=[];
    mindepdXcols=[];
    indepdCrows=[];
    rowstorem=[];
end