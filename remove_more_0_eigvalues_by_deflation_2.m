nullspacesize = size(X,2) - rank(X);
ps=zeros(1,size(X,2));
for p = 1:size(X,2)
t = X; t(:,p)= [];
ps(p) = norm(t*pinv(t)*X(:,p) - X(:,p));
end
depdXcols = find(ps < 1e-10);
C = X(:,depdXcols);
nzrows = find(sum(abs(C))' > 1e-10);
C = C(nzrows,:);

depdCcols = find(diag(rref(C)) <=1e-10)';
mdepdXcols = depdXcols(depdCcols);
indepdCcols = setdiff(1:size(C,2), depdCcols);
mindepdXcols = depdXcols(indepdCcols);
rrefC = rref(C(:, indepdCcols)');
depdCrows = find(sum(abs(rrefC)) == 1);
% Xp1 =  [zeros(size(X,2)-size(C,2),nullspacesize); inv(C(depdCrows,indepdcols))*C(depdCrows,depdCcols); -eye(nullspacesize)];
Xp1 = zeros(size(X,2), nullspacesize);
Xp1(mdepdXcols,:) = -eye(nullspacesize);
Xp1(mindepdXcols,:) = inv(C(depdCrows,indepdCcols))*C(depdCrows,depdCcols);


rows2remove = find( sum(abs(Xp1)') == 1);
Xp2 = X(rows2remove,:);
X = X + Xp1 * Xp2;
Y(mdepdXcols,:)=[];
Y(:,mdepdXcols)=[];
