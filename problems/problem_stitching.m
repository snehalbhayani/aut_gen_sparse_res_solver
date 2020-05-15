function [ eqs, data0, eqs_data ] = problem_stitching2( data0 )

if nargin < 1 || isempty(data0)
    data0 = randi(100,32,1);
end

c2 = data0(17:32);
c1 = data0(1:16);

mm = [ 6     5     4     4     3     3     2     2     2     1     1     1     0     0     0     0;...
     3     3     3     2     3     2     3     2     1     3     2     1     3     2     1     0];
 
eqs = [multipol(c1',mm);multipol(c2',mm)];


if nargout == 3
   xx = create_vars(2+32);
   data = xx(3:end);
   c2 = data(17:32);
   c1 = data(1:16);
   x1 = xx(1);
   x2 = xx(2);
   eqs_data = [sum(c1'.*x1.^mm(1,:).*x2.^mm(2,:));sum(c2'.*x1.^mm(1,:).*x2.^mm(2,:))];
   
end

