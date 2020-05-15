function [eqs,data0,eqs_data] = problem_trisonal_stewenius(data0)

if nargin < 1 || isempty(data0)
    data0 = randi(64,7*4,1);
end

xx = create_vars(4);
Ti = reshape(data0,7,4);
T = Ti*xx;

eq1 = T(5)*T(7)*(2*T(6)+T(7))+T(5)^2*T(7)-T(1);
eq2 = -T(6)*(T(5)+2*T(6)+T(7))-T(2);
eq3 = T(5)*(T(6)+T(7))-T(3);
eq4 = T(7)*(T(5)+T(6))-T(4);

eqs = [eq1; eq2; eq3; eq4];

if nargout == 3
    xx = create_vars(4+7*4);
    data = xx(5:end);
    
    Ti = reshape(data,7,4);
    T = Ti*xx(1:4);
    
    eq1 = T(5)*T(7)*(2*T(6)+T(7))+T(5)^2*T(7)-T(1);
    eq2 = -T(6)*(T(5)+2*T(6)+T(7))-T(2);
    eq3 = T(5)*(T(6)+T(7))-T(3);
    eq4 = T(7)*(T(5)+T(6))-T(4);
    
    eqs_data = [eq1; eq2; eq3; eq4];
    
    totalsyms = 4 + 28;
    noofvars = 4;
    ipparams = strcat('a',num2str(1));
    for i = 2:totalsyms
        if i > noofvars
            ipparams = strjoin({ipparams, ',c', num2str(i-noofvars)}, '');
        else
            ipparams = strjoin({ipparams, ',a', num2str(i)}, '');
        end
    end
    
    fileID = fopen('Eqs_problem_trisonal_stewenius.m','w');
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

