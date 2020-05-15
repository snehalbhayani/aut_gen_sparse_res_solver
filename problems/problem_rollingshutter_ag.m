function [ eqs, data0, eqs_data ] = problem_rollingshutter( data0 )

if nargin < 1 || isempty(data0)
    data0 = randi(50,6*9,1);
 
end

xx = create_vars(5);

rv = xx(1)*data0(1:9)+xx(2)*data0(10:18)+xx(3)*data0(19:27)+xx(4)*data0(28:36)+xx(5)*data0(37:45)+data0(46:54);
r1 = rv(1:3);
r2 = rv(4:6);
r3 = rv(7:9);
c1 = rv([1 4 7]);
c2 = rv([2 5 8]);
c3 = rv([3 6 9]);

eqs = [r1'*r1-r2'*r2;r1'*r1-r3'*r3];
eqs = [eqs;c1'*c1-c2'*c2;c1'*c1-c3'*c3];
eqs = [eqs;r1'*r2;r1'*r3;r2'*r3];
eqs = [eqs;c1'*c2;c1'*c3;c2'*c3];


if nargout == 3
    xx = create_vars(5+6*9);
    data = xx(6:end);


rv = xx(1)*data(1:9)+xx(2)*data(10:18)+xx(3)*data(19:27)+xx(4)*data(28:36)+xx(5)*data(37:45)+data(46:54);
r1 = rv(1:3);
r2 = rv(4:6);
r3 = rv(7:9);
c1 = rv([1 4 7]);
c2 = rv([2 5 8]);
c3 = rv([3 6 9]);

eqs_data = [r1'*r1-r2'*r2;r1'*r1-r3'*r3];
eqs_data = [eqs_data;c1'*c1-c2'*c2;c1'*c1-c3'*c3];
eqs_data = [eqs_data;r1'*r2;r1'*r3;r2'*r3];
eqs_data = [eqs_data;c1'*c2;c1'*c3;c2'*c3];

      totalsyms = 5 + 6*9;
    noofvars = 5;
    ipparams = strcat('a',num2str(1));
    for i = 2:totalsyms
        if i > noofvars
            ipparams = strjoin({ipparams, ',c', num2str(i-noofvars)}, '');
        else
            ipparams = strjoin({ipparams, ',a', num2str(i)}, '');
        end
    end
    
    fileID = fopen('Eqs_problem_rollingshutter.m','w');
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

