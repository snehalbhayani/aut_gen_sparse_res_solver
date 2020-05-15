A=[2, 2; 2, 2]; 
[v,l]=eigs(A,1); 
A=A-l*eye(2,2); 
[v,l]=eigs(A,1); 
[~,p]=max(v); 
%deflation 
A1=A-v/v(p)*A(p,:); 
%removal of row and column 
A1(p,:)=[]; A1(:,p)=[]; 
[v,l]=eigs(A1,1);