function [ eqs, data0, eqs_data ] = problem_toa_55_3( data0 )

if nargin < 1 || isempty(data0)
    data0 = randi(50,71,1);
end

xx = create_vars(6);

vv = [1;xx(1:5)];
lam = xx(6);
H = 0;
b = 0;
count = 0;
for iii = 1:6,
    C = [data0(count+1) data0(count+2) data0(count+3);data0(count+2) data0(count+4) data0(count+5);data0(count+3) data0(count+5) data0(count+6)];
    H = H+vv(iii)*C;
    b = b + data0((7:9)+count)*vv(iii);
    count = count + 9;
end

adjH = adj(H);
detH = det(H);

S = reshape(data0((1:12)+count),3,4);
dd = data0((13:17)+count);
eqs(5,1) = b'*adjH*b-dd(5)*detH;

for iii = 1:4,
    eqs(iii) = (S(:,iii)'*adjH+2*b'*adjH)*S(:,iii)-dd(iii)*detH;
end
eqs(6) = detH-lam;


if nargout > 1
    xx = create_vars(6+71);
    data = xx(7:end);
    
    vv = [1;xx(1:5)];
    lam = xx(6);
    H = 0;
    b = 0;
    count = 0;
    for iii = 1:6,
        C = [data(count+1) data(count+2) data(count+3);data(count+2) data(count+4) data(count+5);data(count+3) data(count+5) data(count+6)];
        H = H+vv(iii)*C;
        b = b + data((7:9)+count)*vv(iii);
        count = count + 9;
    end

    adjH = adj(H);
    detH = det(H);

    S = reshape(data((1:12)+count),3,4);
    dd = data((13:17)+count);
    eqs_data(5,1) = b'*adjH*b-dd(5)*detH;

    for iii = 1:4,
        eqs_data(iii) = (S(:,iii)'*adjH+2*b'*adjH)*S(:,iii)-dd(iii)*detH;
    end
    eqs_data(6) = detH-lam;

end
    
    