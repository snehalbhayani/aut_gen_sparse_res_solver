function [ eqs, data0, eqs_data ] = problem_p4p_fr_iccv17( data0 )
if nargin < 1 || isempty(data0)
    data0 = randi(50,64,1);
end

% Basis giving 28x40 template
% b = [0,2,1,0,0,1,0,1,0,0,0,0;
%      1,0,1,2,1,0,0,0,1,0,0,0;
%      1,0,0,0,1,0,1,0,0,1,0,0;
%      1,0,0,0,0,1,1,0,0,0,1,0]

N = reshape(data0(1:8*4),8,4);
B = reshape(data0(8*4+1:end),4,8);
xx = create_vars(4);
alpha = xx(1:3);
k = xx(4);
P12 = reshape(N*[alpha;1],4,2)';
P3 = (B * [alpha;k*alpha;k;1])';
P = [P12; P3];
x11 = P(1,1); x12 = P(1,2); x13 = P(1,3);
x21 = P(2,1); x22 = P(2,2); x23 = P(2,3);
x31 = P(3,1); x32 = P(3,2); x33 = P(3,3);
eqs(1) = x21*x31+x22*x32+x23*x33;
eqs(2) = x11*x31+x12*x32+x13*x33;
eqs(3) = x11*x21+x12*x22+x13*x23;
eqs(4) = x11^2+x12^2+x13^2-x21^2-x22^2-x23^2;
eqs(5) = x13^2*x32-x21^2*x32-x22^2*x32-x12*x13*x33-x22*x23*x33;
eqs(6) = x12*x13*x32+x22*x23*x32-x12^2*x33+x21^2*x33+x23^2*x33;
eqs(7) = x11*x13*x32+x21*x23*x32-x11*x12*x33-x21*x22*x33;
eqs(8) = x13^2*x31-x22^2*x31+x21*x22*x32-x11*x13*x33;
eqs(9) = x12*x13*x31+x22*x23*x31-x11*x12*x33-x21*x22*x33;

if nargout == 3
    xx = create_vars(4+64);
    eqs_data = problem_p4p_fr_iccv17(xx(5:end));
    
    totalsyms = 4 + 64;
    noofvars = 4;
    ipparams = strcat('a',num2str(1));
    for i = 2:totalsyms
        if i > noofvars
            ipparams = strjoin({ipparams, ',c', num2str(i-noofvars)}, '');
        else
            ipparams = strjoin({ipparams, ',a', num2str(i)}, '');
        end
    end
    
    
    fileID = fopen('Eqs_problem_p4p_fr_iccv17.m','w');
    fprintf(fileID, '%s', strjoin({'function eqs = retrieve_eqs(', ipparams, ') '},''));
    fprintf(fileID, '\n');
    for i = 1:size(eqs,2)
        fprintf(fileID, '%s', strjoin({'eqs(', num2str(i),') = ',char(eqs_data(i), true, [], true, noofvars), ';'},''));
        fprintf(fileID, '\n');
    end
    fprintf(fileID, '\n');
    fprintf(fileID, '%s', 'end');
    fclose(fileID);
    
end