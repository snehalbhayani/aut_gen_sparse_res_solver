function [ eqs, data0, eqs_data ] = problem_optpose3pt( data0 )

if nargin < 1 || isempty(data0)
    data0 = randi(50,9+15,1);
end

U = reshape(data0(1:9),3,3);
C = zeros(3,3,3);
for iii = 1:3,
    v = data0((10:14)+(iii-1)*5);
    Ci = [1 v(1) v(2);v(1) v(3) v(4);v(2) v(4) v(5)];
    C(:,:,iii)=Ci;
end

U = U-repmat(U(:,1),1,3);

xx = create_vars(5);
cth = xx(1);
sth = xx(2);
t = xx(3:5);

R = [1 0 0;0 cth -sth;0 sth cth];
%eqs = ];
eqs = t'*C(:,:,1)*t;
for iii = 2:3,
    Up = R*U(:,iii)+t;
    eqs = [eqs;Up'*C(:,:,iii)*Up];
end
eqs = [eqs;cth^2+sth^2-1];

boffa = [];
for iii = 1:3,
    boffa = [boffa;diff(eqs(iii))];
end
eqs = [eqs;det(boffa(:,3:5))];



if nargout == 3
    xx = create_vars(5+9+15);
    data = xx(6:end);
    

U = reshape(data(1:9),3,3);
Cm(3,3,3)=multipol;
for iii = 1:3,
    v = data((10:14)+(iii-1)*5);
    Ci = [1 v(1) v(2);v(1) v(3) v(4);v(2) v(4) v(5)];
    Cm(:,:,iii)=Ci;
end

U = U-repmat(U(:,1),1,3);

cth = xx(1);
sth = xx(2);
t = xx(3:5);

R = [1 0 0;0 cth -sth;0 sth cth];
%eqs_data = [];
eqs_data = t'*Cm(:,:,1)*t;
for iii = 2:3,
    Up = R*U(:,iii)+t;
    eqs_data = [eqs_data;Up'*Cm(:,:,iii)*Up];
end
eqs_data = [eqs_data;cth^2+sth^2-1];

boffa = [];
for iii = 1:3,
    boffa = [boffa;diff(eqs_data(iii))];
end
eqs_data = [eqs_data;det(boffa(:,3:5))];





end

