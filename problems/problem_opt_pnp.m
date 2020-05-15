function [ eqs, data0, eqs_data ] = problem_opt_pnp( data0 )

if nargin < 1 || isempty(data0)
    data0 = randi(50,20*11,1);
end

M = reshape(data0,20,11);
xx = create_vars(4);
a = xx(1);
b = xx(2);
c = xx(3);
d = xx(4);
alfa = [1 a^2 a*b a*c a*d b^2 b*c b*d c^2 c*d d^2]';
f = alfa'*(M'*M)*alfa;
eqs = diff(f);
eqs = eqs(:);

if nargout == 3
    xx = create_vars(4+20*11);
    data = xx(5:end);
    M = reshape(data,20,11);
    a = xx(1);
    b = xx(2);
    c = xx(3);
    d = xx(4);
    alfa = [1 a^2 a*b a*c a*d b^2 b*c b*d c^2 c*d d^2]';
    f = alfa'*(M'*M)*alfa;
    eqs_data = diff(f);
    eqs_data = eqs_data(:);


end

