function [ eqs, data0, eqs_data ] = problem_opt_pnp_nakanoC( data0 )

if nargin < 1 || isempty(data0)
    M = randi(200,9,9);
    M = M+M';
    data0 = M(:);
    
end

M = reshape(data0,9,9);
xx = create_vars(3);
a = 1;
b = xx(1);
c = xx(2);
d = xx(3);

R = [a^2+b^2-c^2-d^2 2*(b*c-a*d) 2*(b*d+a*c);...
    2*(b*c+a*d) a^2-b^2+c^2-d^2 2*(c*d-a*b);...
    2*(b*d-a*c) 2*(c*d+a*b) a^2-b^2-c^2+d^2];
r = R(:);
matMr = reshape(M*r,3,3);
P = R'*matMr-matMr'*R;
Q = matMr*R'-R*matMr';

eqs = [P(1,2);P(1,3);P(2,3);Q(1,2);Q(1,3);Q(2,3)];


if nargout == 3
    xx = create_vars(3+81);
    eqs_data = problem_opt_pnp_nakanoC(xx(4:end));
    
    totalsyms = 3 + 81;
    noofvars = 3;
    ipparams = strcat('a',num2str(1));
    for i = 2:totalsyms
        if i > noofvars
            ipparams = strjoin({ipparams, ',c', num2str(i-noofvars)}, '');
        else
            ipparams = strjoin({ipparams, ',a', num2str(i)}, '');
        end
    end
    
    fileID = fopen('Eqs_problem_opt_pnp_nakanoC.m','w');
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

