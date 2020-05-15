function [ eqs, data0, eqs_data ] = problem_pose_quiver_orig( data0 )

if nargin < 1 || isempty(data0)
    
    data0 = randi(50,4*9,1);
 
end

xx = create_vars(4);
[M,M2,B]=get_reduced_eqs(data0, xx);
eqs = M*transpose(B) + M2;


if nargout == 3
    xx = create_vars(4+36);
    data = xx(5:end);
    
[M,M2,B]=get_reduced_eqs(data, xx(1:4));
eqs_data = M*transpose(B) + M2;

    totalsyms = 4 + 36;
    noofvars = 4;
    ipparams = strcat('a',num2str(1));
    for i = 2:totalsyms
        if i > noofvars
            ipparams = strjoin({ipparams, ',c', num2str(i-noofvars)}, '');
        else
            ipparams = strjoin({ipparams, ',a', num2str(i)}, '');
        end
    end
    
    fileID = fopen('Eqs_problem_pose_quiver.m','w');
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

