function [ eqs, data0, eqs_data ] = problem_opt_vanishingpoints_sat_eqs( data0 )

if nargin < 1 || isempty(data0)
    data0 = randi(1000,6*3,1);
end

mm1 = data0(1:6);
mm1 = [mm1(1) mm1(2) mm1(3);mm1(2) mm1(4) mm1(5);mm1(3) mm1(5) mm1(6)];
mm2 = data0(7:12);
mm2 = [mm2(1) mm2(2) mm2(3);mm2(2) mm2(4) mm2(5);mm2(3) mm2(5) mm2(6)];
mm3 = data0(13:18);
mm3 = [mm3(1) mm3(2) mm3(3);mm3(2) mm3(4) mm3(5);mm3(3) mm3(5) mm3(6)];

xx = create_vars(4);
a = 1;
b = xx(1);
c = xx(2);
d = xx(3);

s = xx(1:3);
s0 = xx(4);

R = [a^2+b^2-c^2-d^2 2*(b*c-a*d) 2*(b*d+a*c);...
    2*(b*c+a*d) a^2-b^2+c^2-d^2 2*(c*d-a*b);...
    2*(b*d-a*c) 2*(c*d+a*b) a^2-b^2-c^2+d^2];

% Jp = 0;
% for iii = 1:nrp,
%     Jp = Jp + ([1 0 0]*R*mm(:,iii))^2;
% end
% for iii = 1:nrp,
%     Jp = Jp + ([0 1 0]*R*mm(:,nrp+iii))^2;
% end
% for iii = 1:nrp,
%     Jp = Jp + ([0 0 1]*R*mm(:,2*nrp+iii))^2;
% end
% Jp = Jp*0.5;
v1 = [1 0 0]';
v2 = [0 1 0]';
v3 = [0 0 1]';

Jp = (v1'*R*(mm1)*R'*v1)+(v2'*R*(mm2)*R'*v2)+(v3'*R*(mm3)*R'*v3);

dJp = diff(Jp);
eqs = [(1+s'*s)*dJp(1)-4*s(1)*Jp;(1+s'*s)*dJp(2)-4*s(2)*Jp;(1+s'*s)*dJp(3)-4*s(3)*Jp];
eqs = [eqs;s0-(1+s'*s)];


if nargout == 3
    xx = create_vars(4+6*3);
    data = xx(5:end);
    
    
    mm1 = data(1:6);
    mm1 = [mm1(1) mm1(2) mm1(3);mm1(2) mm1(4) mm1(5);mm1(3) mm1(5) mm1(6)];
    mm2 = data(7:12);
    mm2 = [mm2(1) mm2(2) mm2(3);mm2(2) mm2(4) mm2(5);mm2(3) mm2(5) mm2(6)];
    mm3 = data(13:18);
    mm3 = [mm3(1) mm3(2) mm3(3);mm3(2) mm3(4) mm3(5);mm3(3) mm3(5) mm3(6)];
    
    a = 1;
    b = xx(1);
    c = xx(2);
    d = xx(3);
    
    s = xx(1:3);
    s0 = xx(4);
    
    R = [a^2+b^2-c^2-d^2 2*(b*c-a*d) 2*(b*d+a*c);...
        2*(b*c+a*d) a^2-b^2+c^2-d^2 2*(c*d-a*b);...
        2*(b*d-a*c) 2*(c*d+a*b) a^2-b^2-c^2+d^2];
    
    % Jp = 0;
    % for iii = 1:nrp,
    %     Jp = Jp + ([1 0 0]*R*mm(:,iii))^2;
    % end
    % for iii = 1:nrp,
    %     Jp = Jp + ([0 1 0]*R*mm(:,nrp+iii))^2;
    % end
    % for iii = 1:nrp,
    %     Jp = Jp + ([0 0 1]*R*mm(:,2*nrp+iii))^2;
    % end
    % Jp = Jp*0.5;
    v1 = [1 0 0]';
    v2 = [0 1 0]';
    v3 = [0 0 1]';
    
    Jp = (v1'*R*(mm1)*R'*v1)+(v2'*R*(mm2)*R'*v2)+(v3'*R*(mm3)*R'*v3);
    
    dJp = diff(Jp);
    eqs_data = [(1+s'*s)*dJp(1)-4*s(1)*Jp;(1+s'*s)*dJp(2)-4*s(2)*Jp;(1+s'*s)*dJp(3)-4*s(3)*Jp];
    eqs_data = [eqs_data;s0-(1+s'*s)];
    totalsyms = 4 + 18;
    noofvars = 4;
    ipparams = strcat('a',num2str(1));
    for i = 2:totalsyms
        if i > noofvars
            ipparams = strjoin({ipparams, ',c', num2str(i-noofvars)}, '');
        else
            ipparams = strjoin({ipparams, ',a', num2str(i)}, '');
        end
    end
    
    fileID = fopen('Eqs_problem_opt_vanishingpoints_sat_eqs.m','w');
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