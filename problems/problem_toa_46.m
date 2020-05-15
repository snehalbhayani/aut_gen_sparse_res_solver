function [ eqs, data0, eqs_data ] = problem_toa_46( data0 )
if nargin < 1 || isempty(data0)
    zp = 30097;
    r = randi(50,3,6);
    s = randi(50,3,4);
    r(:,1)=[0;0;0];
    d = toa_calc_d_from_xy(r,s);
    d2 = round(d.^2);
    cl = compactionmatrix(6);
    cr = compactionmatrix(4);
    B = cl*d2*cr'/(-2);
    R = B';
    S = eye(3);
    da = d2(1,1);
    db = d2(1,:)*cr';
    dc = cl*d2(:,1);
    H{1} = [1 0 0;0 0 0;0 0 0];
    H{2} = [0 0 0;0 1 0;0 0 0];
    H{3} = [0 0 0;0 0 0;0 0 1];
    H{4} = [0 1 0;1 0 0;0 0 0];
    H{5} = [0 0 1;0 0 0;1 0 0];
    H{6} = [0 0 0;0 0 1;0 1 0];
    b{1} = [1;0;0];
    b{2} = [0;1;0];
    b{3} = [0;0;1];
    % 5 linear equations of type C
    % R(:,k)'*H*R(:,k) -2*b'*R(:,k) = dc(k);
    AA = zeros(5,9);
    bb = zeros(5,1);
    for ii = 1:5;
        for jj = 1:6;
            AA(ii,jj) = R(:,ii)'*H{jj}*R(:,ii);
        end
        for jj = 1:3;
            AA(ii,jj+6) = -2*b{jj}'*R(:,ii);
        end
        bb(ii)=dc(ii);
    end
    
    AA = round(AA);
    bb = round(bb);
    AA = mod(AA,zp);
    bb = mod(bb,zp);
    % Find particular solution zz0 and
    % basis zzb for nullspace
    % so that all solutions can be written
    % zz = zz0 + zzb*xx;
    % for a general vector xx in R^4
    load ik
    [zz0,zzb]=zplineqik(AA,bb,zp,ik);
    xtmp = randi(100,4,5);
    %keyboard;
    zz0 = mod(zz0 + zzb*xtmp(:,5),zp);
    zzb = mod(zzb*xtmp(:,1:4),zp);
    %[xp,xh,sols]=zplineq(AA,bb,zp);
    %     zz0 = AA\bb;
    %     [u,s,v]=svd(AA);
    %     zzb = v(:,6:9);
    data0 = [R(:); zzb(:); zz0(:); db(:); da(:)];
    %data0 = randi(50,64,1);
end

H{1} = [1 0 0;0 0 0;0 0 0];
H{2} = [0 0 0;0 1 0;0 0 0];
H{3} = [0 0 0;0 0 0;0 0 1];
H{4} = [0 1 0;1 0 0;0 0 0];
H{5} = [0 0 1;0 0 0;1 0 0];
H{6} = [0 0 0;0 0 1;0 1 0];
b{1} = [1;0;0];
b{2} = [0;1;0];
b{3} = [0;0;1];


%data0 = randi(20,64,1);
S = eye(3);
R = reshape(data0(1:15),3,5);
zzb = reshape(data0(16:51),9,4);
zz0 = reshape(data0(52:60),9,1);
db = reshape(data0(61:63),3,1);
da = reshape(data0(64),1,1);
xx = create_vars(5);
zz = zz0+zzb*xx(1:4);
HH = zeros(3,3);
for k = 1:6;
    HH = HH + H{k}*zz(k);
end
BB = zeros(3,1);
for k = 1:3;
    BB = BB + b{k}*zz(6+k);
end
eqs(1,1) = S(:,1)'*adj(HH)*S(:,1) + 2*BB'*adj(HH)*S(:,1) - xx(5)*db(1);
eqs(2,1) = S(:,2)'*adj(HH)*S(:,2) + 2*BB'*adj(HH)*S(:,2) - xx(5)*db(2);
eqs(3,1) = S(:,3)'*adj(HH)*S(:,3) + 2*BB'*adj(HH)*S(:,3) - xx(5)*db(3);
eqs(4,1) = BB'*adj(HH)*BB - xx(5)*da;
eqs(5,1) = det(HH)-xx(5);

if nargout > 1
    xx = create_vars(5+64);
    data = xx(6:end);
    
    S = eye(3);
    R = reshape(data(1:15),3,5);
    zzb = reshape(data(16:51),9,4);
    zz0 = reshape(data(52:60),9,1);
    db = reshape(data(61:63),3,1);
    da = reshape(data(64),1,1);
    xx = create_vars(5);
    zz = zz0+zzb*xx(1:4);
    HH = zeros(3,3);
    for k = 1:6;
        HH = HH + H{k}*zz(k);
    end
    BB = zeros(3,1);
    for k = 1:3;
        BB = BB + b{k}*zz(6+k);
    end
    if 0,
        eqs_data(1) = xx(5)*S(:,1)'*adj(HH)*S(:,1) + 2*BB'*adj(HH)*S(:,1) - db(1);
        eqs_data(2) = xx(5)*S(:,2)'*adj(HH)*S(:,2) + 2*BB'*adj(HH)*S(:,2) - db(2);
        eqs_data(3) = xx(5)*S(:,3)'*adj(HH)*S(:,3) + 2*BB'*adj(HH)*S(:,3) - db(3);
        eqs_data(4) = xx(5)*BB'*adj(HH)*BB - da;
        eqs_data(5) = det(HH)-xx(5);
    else
        eqs_data(1) = S(:,1)'*adj(HH)*S(:,1) + 2*BB'*adj(HH)*S(:,1) - xx(5)*db(1);
        eqs_data(2) = S(:,2)'*adj(HH)*S(:,2) + 2*BB'*adj(HH)*S(:,2) - xx(5)*db(2);
        eqs_data(3) = S(:,3)'*adj(HH)*S(:,3) + 2*BB'*adj(HH)*S(:,3) - xx(5)*db(3);
        eqs_data(4) = BB'*adj(HH)*BB - xx(5)*da;
        eqs_data(5) = det(HH)-xx(5);
    end
end
