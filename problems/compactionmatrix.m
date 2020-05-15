function [cc,dd] = compactionmatrix(n);
% cc = compactionmatrix(n);
%

cc = [-ones(n-1,1) eye(n-1)];
dd = [1 zeros(1,n-1);cc];
