function [ eqs, data0, eqs_data ] = problem_optpose2pt_v2( data0 )


if nargin < 1 || isempty(data0)
    data0 = randi(50,6+2*4,1);
end

U = reshape(data0(1:6),3,2);


count = 6;
for iii = 1:2,
    C(:,:,iii) = [data0(1+count) 0 data0(2+count);0 data0(1+count) data0(3+count); data0(2+count) data0(3+count) data0(4+count)];
    count = count+4;
end
    
xx = create_vars(5);
cth = xx(4);
sth = xx(5);
t = xx(1:3);

R = [1 0 0;0 cth -sth;0 sth cth];
eqs = [];
for iii = 1:2,
    Up = R*U(:,iii)+t;
    eqs = [eqs;Up'*C(:,:,iii)*Up];
end
eqs = [eqs;cth^2+sth^2-1];

boffa = [];
for iii = 1:2,
    boffa = [boffa;diff(eqs(iii))];
end

eqs = [eqs;det(boffa(:,[3 4]))];
eqs = [eqs;det(boffa(:,[3 5]))];
eqs = [eqs;det(boffa(:,[4 5]))];

%eqs = [eqs;det([1 0 0 0 0;2*cth 2*sth 0 0 0;boffa])];



if nargout == 3
    xx = create_vars(5+14);
    data = xx(6:end);
    
    
    U = reshape(data(1:6),3,2);
    count = 6;
    for iii = 1:2,
        Cm(:,:,iii) = [data(1+count) 0 data(2+count);0 data(1+count) data(3+count); data(2+count) data(3+count) data(4+count)];
    count = count+4;
   end
 
    cth = xx(4);
    sth = xx(5);
    t = xx(1:3);
    R = [1 0 0;0 cth -sth;0 sth cth];
    eqs_data = [];
    for iii = 1:2,
        Up = R*U(:,iii)+t;
        eqs_data = [eqs_data;Up'*Cm(:,:,iii)*Up];
    end
    eqs_data = [eqs_data;cth^2+sth^2-1];
    boffa = [];
    for iii = 1:2,
        boffa = [boffa;diff(eqs_data(iii))];
    end
    

eqs_data = [eqs_data;det(boffa(:,[3 4]))];
eqs_data = [eqs_data;det(boffa(:,[3 5]))];
eqs_data = [eqs_data;det(boffa(:,[4 5]))];



end

