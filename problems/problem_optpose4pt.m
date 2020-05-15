function [ eqs, data0, eqs_data ] = problem_optpose4pt( data0 )

if nargin < 1 || isempty(data0)
    R0 = [1 0 0;0 3/5 -4/5;0 4/5 3/5];
    t0 = 5*randi(10,3,1);
    e = 1;
    x = [randi(10,2,4);ones(1,4)];
    X0 = (x+[e e e e;0 0 0 0;0 0 0 0]).*repmat(randi(10,1,4),3,1)*5;
    X = round(R0'*X0-repmat(R0'*t0,1,4));
    x = x(1:2,:);
    for (kk = 1:4)        
        C0 = [1 0 -x(1,kk);0 1 -x(2,kk);-x(1,kk) -x(2,kk) x(1:2,kk)'*x(1:2,kk)-e^2];
        C(:,:,kk) = round(C0*randi(10));
    end
    data0 = [X(:);C(:)];
end

U = reshape(data0(1:12),3,4);
C = reshape(data0(13:48),3,3,4);
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
    xx = create_vars(5+12+36);
    data = xx(6:end);
    U = reshape(data(1:12),3,4);
    Cm = reshape(data(13:48),3,3,4);
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

