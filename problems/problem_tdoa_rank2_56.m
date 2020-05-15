function [ eqs, data0, eqs_data ] = problem_tdoa_rank2_56( data0 )

if nargin < 1 || isempty(data0)
    data0 = randi(50,5*6,1);
end


xx = create_vars(6);

m = 5;
n = 6;
rk = 2;
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
        eqs(k,1)=det(c(rrs(k1,:),ccs(k2,:)));
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
            eqs_data(k,1)=det(c(rrs(k1,:),ccs(k2,:)));
        end
    end
end

function [cc,dd] = compactionmatrix(n);
% cc = compactionmatrix(n);
%

cc = [-ones(n-1,1) eye(n-1)];
dd = [1 zeros(1,n-1);cc];

