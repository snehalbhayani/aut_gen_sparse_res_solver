function [ eqs, data0, eqs_data ] = problem_relpose_6p_onefocal( data0 )

if nargin < 1 || isempty(data0)
    data0 = randi(50,3*9,1);
end

xx = create_vars(3);

F0 = reshape(data0(1:9),3,3);
F1 = reshape(data0(10:18),3,3);
F2 = reshape(data0(19:27),3,3);

F = F0 + xx(1)*F1 + xx(2)*F2;
K = [1 0 0;0 1 0;0 0 xx(3)];
E = F*K;


%eqs = (2*(F*K*F')*F*K - sum(diag(F*K*F'))*F*K);
eqs = (2*(E*E')*E - sum(diag(E*E'))*E);
eqs = [eqs(:);det(F)];
for iii = 1:6,
    cc = coeffs(eqs(iii));
    mm = monomials(eqs(iii));
    mm(3,:)=mm(3,:)/2;
    eqs(iii)=multipol(cc,mm);
end
for iii = 7:9,
    cc = coeffs(eqs(iii));
    mm = monomials(eqs(iii));
    mm(3,:)=mm(3,:)-1;
    mm(3,:)=mm(3,:)/2;
    eqs(iii)=multipol(cc,mm);
end


if nargout == 3
    xx = create_vars(3+9*3);
    data = xx(4:end);

    F0 = reshape(data(1:9),3,3);
    F1 = reshape(data(10:18),3,3);
    F2 = reshape(data(19:27),3,3);

    F = F0 + xx(1)*F1 + xx(2)*F2;
    K = [1 0 0;0 1 0;0 0 xx(3)];
    E = F*K;


%eqs = (2*(F*K*F')*F*K - sum(diag(F*K*F'))*F*K);
    eqs_data = (2*(E*E')*E - sum(diag(E*E'))*E);
    eqs_data = [eqs_data(:);det(F)];
    for iii = 1:6,
    cc = coeffs(eqs_data(iii));
    mm = monomials(eqs_data(iii));
    mm(3,:)=mm(3,:)/2;
    eqs_data(iii)=multipol(cc,mm);
end
for iii = 7:9,
    cc = coeffs(eqs_data(iii));
    mm = monomials(eqs_data(iii));
    mm(3,:)=mm(3,:)-1;
    mm(3,:)=mm(3,:)/2;
    eqs_data(iii)=multipol(cc,mm);
end


    
    totalsyms = 3 + 27;
    noofvars = 3;
    ipparams = strcat('a',num2str(1));
    for i = 2:totalsyms
        if i > noofvars
            ipparams = strjoin({ipparams, ',c', num2str(i-noofvars)}, '');
        else
            ipparams = strjoin({ipparams, ',a', num2str(i)}, '');
        end
    end
    
    fileID = fopen('Eqs_problem_relpose_6p_onefocal.m','w');
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

