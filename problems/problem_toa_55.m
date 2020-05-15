function [ eqs, data0, eqs_data ] = problem_toa_55( data0 )
if nargin < 1 || isempty(data0)
    zp = 30097;
    %zp = 1777;
    r = randi(20,3,5);
    s = randi(20,3,5);
    r(:,1)=[0;0;0];
    d = toa_calc_d_from_xy(r,s);
    d2 = round(d.^2);
    cl = compactionmatrix(5);
    cr = compactionmatrix(5);
    B = cl*d2*cr'/(-2);
    R = B(:,1:3)';
    load ik1777
    [zz0,zzb]=zplineqik(mod(R',zp),mod(B(:,4),zp),zp,ik);
    S = [eye(3) zz0];
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
    AA = zeros(4,9);
    bb = zeros(4,1);
    for ii = 1:4;
        for jj = 1:6;
            AA(ii,jj) = R(:,ii)'*H{jj}*R(:,ii);
        end
        for jj = 1:3;
            AA(ii,jj+6) = -2*b{jj}'*R(:,ii);
        end
        bb(ii)=dc(ii);
    end
    
    AA = mod(AA,zp);
    bb = mod(bb,zp);
    AA = round(AA);
    bb = round(bb);
    % Find particular solution zz0 and
    % basis zzb for nullspace
    % so that all solutions can be written
    % zz = zz0 + zzb*xx;
    % for a general vector xx in R^4
    [zz0,zzb]=zplineqik(AA,bb,zp,ik);
    xtmp = randi(20,5,6);
    %keyboard;
    zz0 = mod(zz0 + zzb*xtmp(:,6),zp);
    zzb = mod(zzb*xtmp(:,1:5),zp);
    %[xp,xh,sols]=zplineq(AA,bb,zp);
    %     zz0 = AA\bb;
    %     [u,s,v]=svd(AA);
    %     zzb = v(:,6:9);
    data0 = [R(:); zzb(:); zz0(:); db(:); da(:); S(:,4)];
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
S = [eye(3) data0(72:74)];
R = reshape(data0(1:12),3,4);
zzb = reshape(data0(13:57),9,5);
zz0 = reshape(data0(58:66),9,1);
db = reshape(data0(67:70),4,1);
da = reshape(data0(71),1,1);
xx = create_vars(6);
zz = zz0+zzb*xx(1:5);
HH = zeros(3,3);
for k = 1:6;
    HH = HH + H{k}*zz(k);
end
BB = zeros(3,1);
for k = 1:3;
    BB = BB + b{k}*zz(6+k);
end
eqs(1,1) = S(:,1)'*adj(HH)*S(:,1) + 2*BB'*adj(HH)*S(:,1) - xx(6)*db(1);
eqs(2,1) = S(:,2)'*adj(HH)*S(:,2) + 2*BB'*adj(HH)*S(:,2) - xx(6)*db(2);
eqs(3,1) = S(:,3)'*adj(HH)*S(:,3) + 2*BB'*adj(HH)*S(:,3) - xx(6)*db(3);
eqs(4,1) = S(:,4)'*adj(HH)*S(:,4) + 2*BB'*adj(HH)*S(:,4) - xx(6)*db(4);
eqs(5,1) = BB'*adj(HH)*BB - xx(6)*da;
eqs(6,1) = det(HH)-xx(6);

if nargout > 1
    xx = create_vars(6+74);
    data = xx(7:end);
    
    S = [eye(3) data(72:74)];
    R = reshape(data(1:12),3,4);
    zzb = reshape(data(13:57),9,5);
    zz0 = reshape(data(58:66),9,1);
    db = reshape(data(67:70),4,1);
    da = reshape(data(71),1,1);
    xx = create_vars(6);
    zz = zz0+zzb*xx(1:5);
    HH = zeros(3,3);
    for k = 1:6;
        HH = HH + H{k}*zz(k);
    end
    BB = zeros(3,1);
    for k = 1:3;
        BB = BB + b{k}*zz(6+k);
    end
    if 0,
        eqs_data(1) = xx(6)*S(:,1)'*adj(HH)*S(:,1) + 2*BB'*adj(HH)*S(:,1) - db(1);
        eqs_data(2) = xx(6)*S(:,2)'*adj(HH)*S(:,2) + 2*BB'*adj(HH)*S(:,2) - db(2);
        eqs_data(3) = xx(6)*S(:,3)'*adj(HH)*S(:,3) + 2*BB'*adj(HH)*S(:,3) - db(3);
        eqs_data(4) = xx(6)*S(:,4)'*adj(HH)*S(:,4) + 2*BB'*adj(HH)*S(:,4) - db(4);
        eqs_data(5) = xx(6)*BB'*adj(HH)*BB - da;
        eqs_data(6) = det(HH)-xx(6);
    else
        eqs_data(1) = S(:,1)'*adj(HH)*S(:,1) + 2*BB'*adj(HH)*S(:,1) - xx(6)*db(1);
        eqs_data(2) = S(:,2)'*adj(HH)*S(:,2) + 2*BB'*adj(HH)*S(:,2) - xx(6)*db(2);
        eqs_data(3) = S(:,3)'*adj(HH)*S(:,3) + 2*BB'*adj(HH)*S(:,3) - xx(6)*db(3);
        eqs_data(4) = S(:,4)'*adj(HH)*S(:,4) + 2*BB'*adj(HH)*S(:,4) - xx(6)*db(4);
        eqs_data(5) = BB'*adj(HH)*BB - xx(6)*da;
        eqs_data(6) = det(HH)-xx(6);
    end
end
