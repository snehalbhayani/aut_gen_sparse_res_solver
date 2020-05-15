function [ eqs, data0, eqs_data ] = problem_relpose_6p_rad_1s( data0 )

if nargin < 1 || isempty(data0)
 
    data0 = randi(50,36,1);
 
end

xx = create_vars(3);

e3 = xx(1);
e6 = xx(2);
lam = xx(3);

h1 = data0(1:6);
h2 = data0(7:12);
h3 = data0(13:18);
h4 = data0(19:24);
h5 = data0(25:30);
h6 = data0(31:36);
vv = [lam*e3 lam*e6 e3 e6 lam 1];

e1 = vv*h1;
e2 = vv*h2;
e4 = vv*h3;
e5 = vv*h4;
e7 = vv*h5;
e8 = vv*h6;
e9 = 1;

E = [e1 e4 e7;e2 e5 e8;e3 e6 e9];

eqs = 2*(E*E')*E - sum(diag(E*E'))*E;
eqs = [eqs(:);det(E)];

if nargout == 3
    xx = create_vars(3+36);
    data = xx(4:end);
    
    e3 = xx(1);
    e6 = xx(2);
    lam = xx(3);

    h1 = data(1:6);
    h2 = data(7:12);
    h3 = data(13:18);
    h4 = data(19:24);
    h5 = data(25:30);
    h6 = data(31:36);
    vv = [lam*e3 lam*e6 e3 e6 lam 1];

    e1 = vv*h1;
    e2 = vv*h2;
    e4 = vv*h3;
    e5 = vv*h4;
    e7 = vv*h5;
    e8 = vv*h6;
    e9 = 1;

    E = [e1 e4 e7;e2 e5 e8;e3 e6 e9];

    eqs_data = 2*(E*E')*E - sum(diag(E*E'))*E;
    eqs_data = [eqs_data(:);det(E)];


end

