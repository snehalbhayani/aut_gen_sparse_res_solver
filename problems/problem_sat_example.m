function [ eqs, data0, eqs_data ] = problem_sat_example( data0 )
if nargin < 1 || isempty(data0)
    data0 = randi(50,30,1);
end


xx = create_vars(3);
eqs1 = [xx(1); xx(2:3)'*xx(2:3)-1];
mm = [xx(1);xx(2);xx(3);1];
A = mm*mm';
mm = A([1 5 6 9 10 11 13 14 15 16])';
eqs2 = reshape(data0,3,10)*mm;
eqs = eqs1*eqs2';
eqs = eqs(:);

if nargout > 1
    xx = create_vars(3+30);
    data = xx(4:end);
    eqs1 = [xx(1); xx(2:3)'*xx(2:3)-1];
    mm = [xx(1);xx(2);xx(3);1];
    A = mm*mm';
    mm = A([1 5 6 9 10 11 13 14 15 16])';
    eqs2 = reshape(data,3,10)*mm;
    eqs_data = eqs1*eqs2';
    eqs_data = eqs_data(:);
end
    


end

