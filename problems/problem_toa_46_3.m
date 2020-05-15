function [ eqs, data0, eqs_data ] = problem_toa_46_3( data0 )

if nargin < 1 || isempty(data0)
    data0 = randi(50,58,1);
end

xx = create_vars(5);

vv = [1;xx(1:4)];
lam = xx(5);
H = 0;
b = 0;
count = 0;
for iii = 1:5,
    C = [data0(count+1) data0(count+2) data0(count+3);data0(count+2) data0(count+4) data0(count+5);data0(count+3) data0(count+5) data0(count+6)];
    H = H+vv(iii)*C;
    b = b + data0((7:9)+count)*vv(iii);
    count = count + 9;
end

adjH = adj(H);
detH = det(H);

S = reshape(data0((1:9)+count),3,3);
dd = data0((10:13)+count);
eqs(4,1) = b'*adjH*b-dd(4)*detH;

for iii = 1:3,
    eqs(iii) = (S(:,iii)'*adjH+2*b'*adjH)*S(:,iii)-dd(iii)*detH;
end
eqs(5) = detH-lam;


if nargout > 1
    xx = create_vars(5+58);
    data = xx(6:end);
    
    vv = [1;xx(1:4)];
    lam = xx(5);
    H = 0;
    b = 0;
    count = 0;
    for iii = 1:5,
        C = [data(count+1) data(count+2) data(count+3);data(count+2) data(count+4) data(count+5);data(count+3) data(count+5) data(count+6)];
        H = H+vv(iii)*C;
        b = b + data((7:9)+count)*vv(iii);
        count = count + 9;
    end

    adjH = adj(H);
    detH = det(H);

    S = reshape(data((1:9)+count),3,3);
    dd = data((10:13)+count);
    eqs_data(4,1) = b'*adjH*b-dd(4)*detH;

    for iii = 1:3,
        eqs_data(iii) = (S(:,iii)'*adjH+2*b'*adjH)*S(:,iii)-dd(iii)*detH;
    end
    eqs_data(5) = detH-lam;

end
    
    