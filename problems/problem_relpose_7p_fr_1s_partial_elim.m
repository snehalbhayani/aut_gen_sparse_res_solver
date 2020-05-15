function [ eqs, data0, eqs_data ] = problem_relpose_7p_fr_1s_partial_elim( data0 )
int_data = 0;
if nargin < 1
    int_data = 1;
    xx = create_vars(8+1);

    F = reshape([xx(1:8);1],3,3);
    k = xx(9);

    F = [F F(:,3)*k];

    x1 = randi(100,2,7);
    x2 = randi(100,2,7);
    d1 = sum(x1.^2);

    eqs = diag([x2;ones(1,7)]'*F*[x1;ones(1,7);d1]);

    [cc,mm] = polynomials2matrix(eqs);

    ind = [1 3 4 5 6 7 8];
    ind0 = setdiff(1:12,ind);

    A = cc(:,ind);
    B = cc(:,ind0);

    G = -zp_matrix_inverse(A,30097)*B;
    G = mod(G,30097);
    data0 = G';
    data0 = data0(:);
end

G = reshape(data0,5,7)';

xx = (create_vars(3));

f13 = xx(1);
f23 = xx(2);
lambda = xx(3);

vv = [f23*lambda f13 f23 lambda 1].';

lf13 = G(1,:)*vv;
f11 = G(2,:)*vv;
f21 = G(3,:)*vv;
f31 = G(4,:)*vv;
f12 = G(5,:)*vv;
f22 = G(6,:)*vv;
f32 = G(7,:)*vv;
f33 = 1;

F = [f11 f12 f13; f21 f22 f23; f31 f32 f33];


eqs = [lf13 - lambda*f13;f13*f22*f31-f12*f23*f31-f13*f21*f32+f11*f23*f32+f12*f21*f33-f11*f22*f33;...
      f11*f13*f23*f31+f21*f23*f23*f31+f12*f13*f23*f32+f22*f23*f23*f32-f11*f13*f21*f33-f12*f13*f22*f33-f21*f21*f23*f33-f22*f22*f23*f33+f23*f31*f31*f33+f23*f32*f32*f33-f21*f31*f33*f33-f22*f32*f33*f33;...
      f11*f13*f13*f31+f13*f21*f23*f31+f12*f13*f13*f32+f13*f22*f23*f32-f11*f11*f13*f33-f12*f12*f13*f33-f11*f21*f23*f33-f12*f22*f23*f33+f13*f31*f31*f33+f13*f32*f32*f33-f11*f31*f33*f33-f12*f32*f33*f33;...
      f11*f13*f13*f21+f12*f13*f13*f22-f11*f11*f13*f23-f12*f12*f13*f23+f13*f21*f21*f23+f13*f22*f22*f23-f11*f21*f23*f23-f12*f22*f23*f23+f13*f21*f31*f33-f11*f23*f31*f33+f13*f22*f32*f33-f12*f23*f32*f33];
% eqs = simplify(eqs);
if int_data
    eqs = zp_reduce(eqs,30097);
end
eqs = multipol(eqs);

if nargout == 3
    xx = create_vars(3+7*5);
    data = xx(4:end);
    G = reshape(data,5,7)';

    f13 = xx(1);
    f23 = xx(2);
    lambda = xx(3);

    vv = [f23*lambda f13 f23 lambda 1].';

    lf13 = G(1,:)*vv;
    f11 = G(2,:)*vv;
    f21 = G(3,:)*vv;
    f31 = G(4,:)*vv;
    f12 = G(5,:)*vv;
    f22 = G(6,:)*vv;
    f32 = G(7,:)*vv;
    f33 = 1;

    F = [f11 f12 f13; f21 f22 f23; f31 f32 f33];

    eqs_data = [lf13 - lambda*f13;f13*f22*f31-f12*f23*f31-f13*f21*f32+f11*f23*f32+f12*f21*f33-f11*f22*f33;...
      f11*f13*f23*f31+f21*f23*f23*f31+f12*f13*f23*f32+f22*f23*f23*f32-f11*f13*f21*f33-f12*f13*f22*f33-f21*f21*f23*f33-f22*f22*f23*f33+f23*f31*f31*f33+f23*f32*f32*f33-f21*f31*f33*f33-f22*f32*f33*f33;...
      f11*f13*f13*f31+f13*f21*f23*f31+f12*f13*f13*f32+f13*f22*f23*f32-f11*f11*f13*f33-f12*f12*f13*f33-f11*f21*f23*f33-f12*f22*f23*f33+f13*f31*f31*f33+f13*f32*f32*f33-f11*f31*f33*f33-f12*f32*f33*f33;...
      f11*f13*f13*f21+f12*f13*f13*f22-f11*f11*f13*f23-f12*f12*f13*f23+f13*f21*f21*f23+f13*f22*f22*f23-f11*f21*f23*f23-f12*f22*f23*f23+f13*f21*f31*f33-f11*f23*f31*f33+f13*f22*f32*f33-f12*f23*f32*f33];
      
          
    totalsyms = 3 + 35;
    noofvars = 3;
    ipparams = strcat('a',num2str(1));
    for i = 2:totalsyms
        if i > noofvars
            ipparams = strjoin({ipparams, ',c', num2str(i-noofvars)}, '');
        else
            ipparams = strjoin({ipparams, ',a', num2str(i)}, '');
        end
    end
    
    fileID = fopen('Eqs_problem_relpose_7p_fr_1s_partial_elim.m','w');
    fprintf(fileID, '%s', strjoin({'function eqs = retrieve_eqs(', ipparams, ') '},''));
    fprintf(fileID, '\n');
    for i = 1:size(eqs,1)
        fprintf(fileID, '%s', strjoin({'eqs(', num2str(i),') = ',char(eqs_data(i,1), true, [], true, noofvars), ';'},''));
        fprintf(fileID, '\n');
    end
    fprintf(fileID, '\n');
    fprintf(fileID, '%s', 'end');
    fclose(fileID);
end



end

