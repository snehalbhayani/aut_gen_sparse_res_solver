function [ eqs, data0, eqs_data ] = problem_optpose4pt_v2( data0 )

if nargin < 1 || isempty(data0)
    data0 = randi(50,12+4*4,1);
end

U = reshape(data0(1:12),3,4);


count = 12;
for iii = 1:4,
    C(:,:,iii) = [data0(1+count) 0 data0(2+count);0 data0(1+count) data0(3+count); data0(2+count) data0(3+count) data0(4+count)];
    count = count+4;
end
    
xx = create_vars(5);
cth = xx(4);
sth = xx(5);
t = xx(1:3);

R = [1 0 0;0 cth -sth;0 sth cth];
eqs = [];
for iii = 1:4,
    Up = R*U(:,iii)+t;
    eqs = [eqs;Up'*C(:,:,iii)*Up];
end
eqs = [eqs;cth^2+sth^2-1];


if nargout == 3
    xx = create_vars(5+28);
    data = xx(6:end);
    
    
    U = reshape(data(1:12),3,4);
    count = 12;
    for iii = 1:4,
        Cm(:,:,iii) = [data(1+count) 0 data(2+count);0 data(1+count) data(3+count); data(2+count) data(3+count) data(4+count)];
    count = count+4;
   end
 
    cth = xx(4);
    sth = xx(5);
    t = xx(1:3);
    R = [1 0 0;0 cth -sth;0 sth cth];
    eqs_data = [];
    for iii = 1:4,
        Up = R*U(:,iii)+t;
        eqs_data = [eqs_data;Up'*Cm(:,:,iii)*Up];
    end
    eqs_data = [eqs_data;cth^2+sth^2-1];

end

