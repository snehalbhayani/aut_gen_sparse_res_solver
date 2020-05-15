function [ eqs, data0, eqs_data ] = problem_unsynch_relpose( data0 )

if nargin < 1 || isempty(data0)
    data0 = randi(50,15*7,1);
end

xx = create_vars(7);
bb = xx(7);
nn = reshape(data0,15,7);
aa = xx(1:6);
w = nn*[aa;1];

eqs(7,1) = det([w(1) w(2) w(3);w(4) w(5) w(6) ;w(7) w(8) w(9)]);
for iii = 1:6,  
    eqs(iii) = w(iii)*bb-w(iii+9);
end

if nargout == 3
    xx = create_vars(7+15*7);
    data = xx(8:end);
    bb = xx(7);
    nn = reshape(data,15,7);
    aa = xx(1:6);
    w = nn*[aa;1];

    eqs_data(7,1) = det([w(1) w(2) w(3);w(4) w(5) w(6) ;w(7) w(8) w(9)]);
    for iii = 1:6,  
        eqs_data(iii) = w(iii)*bb-w(iii+9);
    end
   
    
end

