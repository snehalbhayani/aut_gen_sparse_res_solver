function [ eqs, data0, eqs_data ] = problem_tdoa_rank3_95( data0 )

m = 9;
n = 5;
if nargin < 1 || isempty(data0)
    data0 = randi(50,m*n,1);
end


xx = create_vars(n);

rk = 3;
d = reshape(data0(1:(m*n)),m,n);
[cc1,dd1] = compactionmatrix(m);
[cc2,dd2] = compactionmatrix(n);
tmp1a = (d).^2;
tmp1b = -2*(d.*repmat(xx',m,1));
tmp2a = cc1*tmp1a*cc2'; % This is a consant matrix
tmp2b = cc1*tmp1b*cc2'; % This depends on o
c = tmp2a+tmp2b;

rrs = nchoosek(1:(m-1),rk+1);
ccs = nchoosek(1:(n-1),rk+1);
%rrs = rrs(1:2,:);
%ccs = ccs(1:3,:);
k = 0;
clear eqs;
for k1 = 1:size(rrs,1);
    for k2 = 1:size(ccs,1);
        k = k+1;
        tmp = det(c(rrs(k1,:),ccs(k2,:)));
        ccc = coeffs(tmp);
        mmm = monomials(tmp);
        oki = find(max(mmm) <= 1);
        tmp = multipol(ccc(oki),mmm(:,oki));
        eqs(k,1)=tmp;
    end
end

if nargout == 3
    xx = create_vars(n+m*n);
    oo = xx(1:n);
    data = xx((n+1):end);
    
    
    d = reshape(data(1:(m*n)),m,n);
    [cc1,dd1] = compactionmatrix(m);
    [cc2,dd2] = compactionmatrix(n);
    tmp1a = (d).^2;
    tmp1b = -2*(d.*repmat(xx',m,1));
    tmp2a = cc1*tmp1a*cc2'; % This is a consant matrix
    tmp2b = cc1*tmp1b*cc2'; % This depends on o
    c = tmp2a+tmp2b;
    
    k = 0;
    clear eqs_data;
    for k1 = 1:size(rrs,1);
        for k2 = 1:size(ccs,1);
            k = k+1;
            tmp = det(c(rrs(k1,:),ccs(k2,:)));
            ccc = coeffs(tmp);
            mmm = monomials(tmp);
            oki = find(max(mmm(1:5,:)) <= 1);
            tmp = multipol(ccc(oki),mmm(:,oki));
            eqs_data(k,1)=tmp;
        end
    end
end

function [cc,dd] = compactionmatrix(n);
% cc = compactionmatrix(n);
%

cc = [-ones(n-1,1) eye(n-1)];
dd = [1 zeros(1,n-1);cc];



