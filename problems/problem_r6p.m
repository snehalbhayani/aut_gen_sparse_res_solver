function [ eqs, data0, eqs_data ] = problem_r6p( data0 )


skew = @(v) cross(v(:)*[1 1 1],eye(3));
r0 = 0;

if nargin < 1 || isempty(data0)
    % generate integer instance (somewhat jobbigt)
    data0 = generate_integer_instance();
end



x = [reshape(data0(1:6),2,3);ones(1,3)];
X = reshape(data0(7:15),3,3);
A = reshape(data0(16:end),6,16);

xx = create_vars(3+3);
v = xx(1:3);
w = xx(4:6);

mm = [xx(1)*xx(4) xx(2)*xx(4) xx(3)*xx(4) xx(1)*xx(5) xx(2)*xx(5) xx(3)*xx(5) xx(1)*xx(6) xx(2)*xx(6) xx(3)*xx(6) xx' 1]';
C = A(1:3,:)*mm;
t = A(4:6,:)*mm;

eqs = [];
for k = 1:3
    dr = x(1,k)-r0;
    sx = skew(x(:,k));
    rhs = (eye(3)+dr*skew(w))*(eye(3)+skew(v))*X(:,k)+C+dr*t;
    eqs = [eqs; sx(1:2,:)*rhs];
end


if nargout == 3
    xx = create_vars(3+3+111);
    data = xx(7:end);
    x = [reshape(data(1:6),2,3);ones(1,3)];
    X = reshape(data(7:15),3,3);
    A = reshape(data(16:end),6,16);

    v = xx(1:3);
    w = xx(4:6);

    mm = [xx(1)*xx(4) xx(2)*xx(4) xx(3)*xx(4) xx(1)*xx(5) xx(2)*xx(5) xx(3)*xx(5) xx(1)*xx(6) xx(2)*xx(6) xx(3)*xx(6) xx(1:6)' 1]';
    C = A(1:3,:)*mm;
    t = A(4:6,:)*mm;

    eqs_data = [];
    for k = 1:3
        dr = x(1,k)-r0;
        sx = skew(x(:,k));
        rhs = (eye(3)+dr*skew(w))*(eye(3)+skew(v))*X(:,k)+C+dr*t;
        eqs_data = [eqs_data; sx(1:2,:)*rhs];
    end
end

end

function data0 = generate_integer_instance()
    skew = @(v) cross(v(:)*[1 1 1],eye(3));
    r0 = 0;

    p = 30097;
    % using 3 points we solve for C and t
    x = [randi(50,2,3); 1 1 1];
    X = randi(50,3,3);
    
    xx = create_vars(3+3+3+3);
    v = xx(1:3);
    w = xx(4:6);
    C = xx(7:9);
    t = xx(10:12);

    eqs = [];
    for k = 1:3
        dr = x(1,k)-r0;
        sx = skew(x(:,k));
        rhs = (eye(3)+dr*skew(w))*(eye(3)+skew(v))*X(:,k)+C+dr*t;
        eqs = [eqs; sx(1:2,:)*rhs];
    end
    
    [CC,~] = polynomials2matrix(eqs);

    ind_uw = [1:15 22];
    ind_ct = [16:21];

    if det(CC(:,ind_ct)) == 0
        % degenerate, try again
        data0 = generate_integer_instance();
        return;
    end
    
    Cinv = zp_matrix_inverse(CC(:,ind_ct),p);
    A = mod(Cinv*CC(1:6,ind_uw),p);
    % [C;t] = A*mm;
    
    % data is 3 new random points and the mapping A
    data0 = [randi(50,6+9,1); A(:)];
end
