function [ eqs, data0, eqs_data ] = problem_gp4p_scale( data0 )

if nargin < 1 || isempty(data0)
    data0 = randi(50,13*6,1);
end

b = reshape(data0,13,6);

xx = create_vars(5);

bb = b*[xx;1];

R = reshape(bb(1:9),3,3)';

len_c = sum(R.^2,1);
len_r = sum(R.^2,2);

eqs = [len_r(1)-len_r(2);
       len_r(1)-len_r(3);
       len_c(1)-len_r(2);
       len_c(1)-len_c(3);
       R(1,:)*R(2,:)';
       R(1,:)*R(3,:)';
       R(2,:)*R(3,:)';
       R(:,1)'*R(:,2);
       R(:,1)'*R(:,3);
       R(:,2)'*R(:,3)];


if nargout == 3


    xx = create_vars(5+13*6);
    
    data = xx(6:end);
    b = reshape(data,13,6);

    bb = b*[xx(1:5);1];

    R = reshape(bb(1:9),3,3)';

    len_c = sum(R.^2,1);
    len_r = sum(R.^2,2);

    eqs_data = [len_r(1)-len_r(2);
           len_r(1)-len_r(3);
           len_c(1)-len_r(2);
           len_c(1)-len_c(3);
           R(1,:)*R(2,:)';
           R(1,:)*R(3,:)';
           R(2,:)*R(3,:)';
           R(:,1)'*R(:,2);
           R(:,1)'*R(:,3);
           R(:,2)'*R(:,3)];

       
      totalsyms = 5 + 13*6;
    noofvars = 5;
    ipparams = strcat('a',num2str(1));
    for i = 2:totalsyms
        if i > noofvars
            ipparams = strjoin({ipparams, ',c', num2str(i-noofvars)}, '');
        else
            ipparams = strjoin({ipparams, ',a', num2str(i)}, '');
        end
    end
    
    fileID = fopen('Eqs_problem_gp4p_scale.m','w');
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

