function [q_gt, t_gt, r_gt, datas_p5p_refractive] = generate_3D_p5p_refractive_scene(iters)
clc;
%% A synthetic scene setup for 2D case.
datas_p5p_refractive = [];
syms vx vy vz;
v_sym = [vx;vy;vz];
thresh = 1e-10;
%% Generating 500 differnt data points

%% Ground truth values for the camera pose
% t_gt = [randi(97,2,1);0];
q_gt = [1;randn(3,1)];
R = quat2rot(q_gt);
alpha = norm(R); 
R = R/norm(R);

t_gt = cross(R(:,3), randn(3,1));
% t_gt = randn(3,1);
t = t_gt;
C = t;

nc = [0;0;1];
n = R * nc;

d = randn();
r = 0.3;
r_gt = r;
impoints = 2 * (rand(iters, 15, 2) - 0.5);

for iter = 1:iters
datasm = [];

impointsbatch = reshape(impoints(iter,:,:),15,2);
for i = 1:15
% while size(datasm,1) < 15 
    %%  Measurements
    
    yw = R*[transpose(impointsbatch(i,:));-(d+transpose(n)*t)] + t;
    
%     u = randn(3,1); u = u/norm(u); % Random image data point representative in world coordinate system
    u = yw - C; u = u/norm(u);
    tht1 = abs(acos(dot(n, u)) - pi);
    tht2 = asin(sin(tht1)*r);
    P = yw;
    %     P = C - u * (transpose(nc) * transpose(R) * C + d)/(transpose(nc) * transpose(R) * u);

%     while (P(1)-C(1))/u(1) < 0 || abs(sin(tht1)*r) > 1
%         u = randn(3,1); u = u/norm(u); % Random image data point representative in world coordinate system
%         tht1 = abs(acos(dot(n, u)) - pi);
%         tht2 = asin(sin(tht1)*r);
%         P = C - u * (transpose(nc) * transpose(R) * C + d)/(transpose(nc) * transpose(R) * u);
%     end
    
    %% Bulding up the eqs to extract the actual rays in world coordinates
    
    eq1 = (transpose(n) * v_sym)^2 - (cos(tht2))^2;
    eq2 = transpose(cross(u,n)) * v_sym;
    eq3 = vx^2 + vy^2 + vz^2 - 1;
    [x,y,z] = solve([eq1;eq2;eq3], [vx;vy;vz]);
    for i=1:length(x)
        v = eval([x(i);y(i);z(i)]);
        X = P - rand() * v;
        if (transpose(n)*C+d)*(transpose(n)*X+d) < 0 && norm(r * cross(u,n) - cross(v,n)) < thresh
            break;
        else
            X = [];
            v = [];
        end
    end
    
    if ~isempty(X) && isreal(X)
        datasm = [datasm;[u, transpose(R)*X * alpha]];
    else
        fprintf(" failed to setup \n");
    end
    
    if size(datasm,1) == 15
        break;
    end    
end
datas_p5p_refractive = [datas_p5p_refractive; [transpose(datasm(:,1)), transpose(datasm(:,2))]];
% %%reverse coding now for the actual equations
% clearvars -except q_gt t_gt r_gt datas d Ps;
% syms tx ty tz q1 q2 q3
% t = [tx;ty;tz];
% q = [1;q1;q2;q3];
% R = quat2rot(q);
% C = t;
% nc = [0;0;1];
% n = R * nc;
% r = r_gt;
% eqs = [];
u = reshape(datas_p5p_refractive(end, 1:15),3,5);
X = reshape(datas_p5p_refractive(end, 16:end),3,5);

transpose(cross(alpha*transpose(R)*u,X-alpha*transpose(R)*C))*[0;0;1]
    
%     for k=1:5
%     u = U(:,k);    
%     eqs = [eqs; transpose(n) * cross(u, X(:,k)-C)];
%     end
disp(iter);

end
gt_p5p = {q_gt, t_gt, r_gt};
save("gt_p5p", "gt_p5p");
save("datas_p5p_refractive","datas_p5p_refractive");

end
