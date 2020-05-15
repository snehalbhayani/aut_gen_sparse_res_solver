function [ eqs, data0, eqs_data ] = problem_relpose_4pt( data0 )

if nargin < 1 || isempty(data0)
 
    data0 = [randi(50,16,1);3/5;4/5];
 
end

x1 = [reshape(data0(1:8),2,4);ones(1,4)];
x2 = [reshape(data0(9:16),2,4);ones(1,4)];
cth = data0(17);
sth = data0(18);
xx = create_vars(5);

tx = [0 -1 xx(2);1 0 -xx(1);-xx(2) xx(1) 0];
r = xx(3:5);
E = 5*tx*(cth*eye(3) +(1-cth)*(r*r')+sth*[0 -r(3) r(2);r(3) 0 -r(1);-r(2) r(1) 0]);

eqs(5,1) = multipol;
for iii = 1:4,
    eqs(iii) = x1(:,iii)'*E*x2(:,iii);
end

eqs(5)=r'*r-1;




if nargout == 3
    xx = create_vars(5+18);
    data = xx(6:end);
    

x1 = [reshape(data(1:8),2,4);ones(1,4)];
x2 = [reshape(data(9:16),2,4);ones(1,4)];
cth = data(17);
sth = data(18);


tx = [0 -1 xx(2);1 0 -xx(1);-xx(2) xx(1) 0];
r = xx(3:5);
E = 5*tx*(cth*eye(3) +(1-cth)*(r*r')+sth*[0 -r(3) r(2);r(3) 0 -r(1);-r(2) r(1) 0]);

eqs_data(5,1) = multipol;
for iii = 1:4,
    eqs_data(iii) = x1(:,iii)'*E*x2(:,iii);
end

eqs_data(5)=r'*r-1;


    totalsyms = 5 + 18;
    noofvars = 5;
    ipparams = strcat('a',num2str(1));
    for i = 2:totalsyms
        if i > noofvars
            ipparams = strjoin({ipparams, ',c', num2str(i-noofvars)}, '');
        else
            ipparams = strjoin({ipparams, ',a', num2str(i)}, '');
        end
    end
    
    fileID = fopen('Eqs_problem_relpose_4pt.m','w');
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

