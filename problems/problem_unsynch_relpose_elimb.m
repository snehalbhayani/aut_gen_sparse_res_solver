function [ eqs, data0, eqs_data ] = problem_unsynch_relpose_elimb( data0 )

if nargin < 1 || isempty(data0)
    data0 = randi(50,15*7,1);
end

xx = create_vars(6);
nn = reshape(data0,15,7);
aa = xx(1:6);
w = nn*[aa;1];

boffa(7,2) = det([w(1) w(2) w(3);w(4) w(5) w(6) ;w(7) w(8) w(9)]);
for iii = 1:6,  
    boffa(iii,1) = w(iii);
    boffa(iii,2) = -w(iii+9);
end
ids = nchoosek(1:7,2);
N = size(ids,1);
eqs(N+1,1)=det([w(1) w(2) w(3);w(4) w(5) w(6) ;w(7) w(8) w(9)]);
for iii = 1:N,
    eqs(iii,1) = det(boffa(ids(iii,:),:));
end
%eqs(end+1) = xx(7)-w(1);

if nargout == 3
    xx = create_vars(7+15*7);
    data = xx(8:end);
    
    nn = reshape(data,15,7);
    aa = xx(1:6);
    w = nn*[aa;1];

    boffa(6,2) = multipol;
    for iii = 1:6,  
        boffa(iii,1) = w(iii);
        boffa(iii,2) = -w(iii+9);
    end
    ids = nchoosek(1:6,2);
    N = size(ids,1);
    eqs_data(N+1,1)=det([w(1) w(2) w(3);w(4) w(5) w(6) ;w(7) w(8) w(9)]);
    for iii = 1:N,
        eqs_data(iii) = det(boffa(ids(iii,:),:));
    end

    eqs_data(end+1) = xx(7)-w(1);

   
    
end

