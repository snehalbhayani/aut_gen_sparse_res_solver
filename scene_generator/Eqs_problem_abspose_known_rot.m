function eqs = retrieve_eqs(a1,a2,a3,a4,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,varargin)
thresh = 1e-10;

% d = 0.5;
% r = 0.5;
% syms c13 c14
%%
% ------------ TRUE VALUE TESTER
a1_t = 0.1;
a2_t = 0.3;
a3_t = 0.1;
a4_t = 0.5;
c13_t = 1.2;
c14_t = 1.5;
c15_t = 0.5;
c16_t = 0.6;
% d = c15_t;
% r= c16_t;
% ------------ TRUE VALUE TESTER
d = c15;
r = c16;
% c13 = 1.2;
% c14 = 1.3;
C = [a1;a2;a3];
% a4 = tan(qy/2);
% c13 = tan(qx/2)
cx = (1-c13 ^2);
sx = 2* c13;
cz = (1-c14 ^2);
sz = 2* c14;

Rx = [[(1+c13^2),0,0];[0,cx,sx];[0,-sx,cx]];
Rz = [[cz,sz,0];[-sz,cz,0];[0,0,(1+c14^2)]];

c = (1-a4 ^2);
s = 2* a4;
R =  Rx * Rz * [[c,0,-s];[0,(1+a4^2),0];[s,0,c]];
n = simplify(R * [0;0;1]);
u1 = [c1;c2;c3];
u2 = [c4;c5;c6];

X1 = [c7;c8;c9];
X2 = [c10;c11;c12];
n_t = eval(subs(n, [a4;c13;c14], [a4_t; c13_t; c14_t]));
C_t = eval(subs(C, [a1;a2;a3;a4], [a1_t;a2_t;a3_t;a4_t]));
%% POint 1
u=u1;
ntu = (transpose(n)*u);
P = C - u * (transpose(n)*C + d)/ntu;
% -------------------------------------------------------------------------
% X1_t = X1;
% -------------------------------------------------------------------------    
v = (X1 - P);v=ntu*v;
% v = simplify((X1_t - P));
%%
part11 = r^2 * (transpose(v)*v) * cross(u,n) .^2;
part12 = (cross(v,n) .^ 2);
eqs(1:3) = (part11 - part12);
% eval(subs(part11 - part12, [a1;a2;a3;a4], [a1_t;a2_t;a3_t;a4_t]))
% ##############################################################################################
%% POint 2
u=u2;
ntu = (transpose(n)*u);
P = C - u * (transpose(n)*C + d)/ntu;
% -------------------------------------------------------------------------
v = (X2 - P);v=ntu*v;
% v = simplify((X1_t - P));
%%
part21 = r^2 * (transpose(v)*v) * cross(u,n) .^2;
part22 = (cross(v,n) .^ 2);
eqs(4:6) = (part21 - part22);
% ############################################################################################################


%% Making of the constraint equations
eqs(7) = transpose(n) * cross(u1, X1-C) ;
eqs(8) = transpose(n) * cross(u2, X2-C);
[N, D] = numden(eqs);
eqs = N;
%% GEnerating GT 
% ---------------------------------------------
if ~isempty(varargin)
    datas=[];
    for t = 1:varargin{1}
        
        [u1_t, X1_t, P1_t] = find_v(n_t, c16_t, C_t, c13_t, c14_t, thresh, c15_t);
        [u2_t, X2_t, P2_t] = find_v(n_t, c16_t, C_t, c13_t, c14_t, thresh, c15_t);
        if ~isempty(X1_t) && ~isempty(X2_t)
            datas = [datas; transpose([u1_t;u2_t;X1_t;X2_t;c13_t;c14_t;c15_t;c16_t])];
        end
    end
    save("problem_abspose_known_rot_refractive","datas");
end
end