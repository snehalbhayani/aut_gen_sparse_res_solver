clc;
clear;
clearvars;

syms l1 l2 l3 l4;
vars = [l1,l2,l3];
hiddenVar = l3;
theta1 = rand();
theta2 = rand();
theta3 = rand();
R = [1 0 0; 0 cos(theta1) -sin(theta1); 0 sin(theta1) cos(theta1)] * [cos(theta2) 0 sin(theta2); 0 1 0; -sin(theta2) 0 cos(theta2)] * [cos(theta3) -sin(theta3) 0; sin(theta3) cos(theta3) 0; 0 0 1];
t = [randi(30);randi(10);1];
focalLen = randi(20);
K_true = diag([focalLen, focalLen, 1]);
tx = [0, -t(3), t(2); t(3),0,-t(1);-t(2),t(1),0];
F_true = inv(K_true) * tx * R * inv(K_true);
% E_true = E_true/norm(E_true,'fro');
Fv_true = reshape(F_true,9,1);
P_true = K_true * [R t];
% The world coordinates
Z_true = [randn(3,6); ones(1,6)];
Y_true = (P_true * Z_true);
X_true = (K_true *[eye(3), zeros(3,1)] * Z_true);
% Hence we will have Y' * E_true * X = 0
Y = bsxfun(@rdivide, Y_true,Y_true(3,:));
X = bsxfun(@rdivide, X_true,X_true(3,:));
% X = Z_true(1:3,:) ./ Z_true(3,:);

% This is the random part. We assume 5 - random point correspondences in
% image coordinate system.
qs = zeros(6,9);
% With the point correspondences known, we estimate a parametric form of
% essential matrix as a linear combination of the null space of a function
% of the pair of corresponding points, qs = [x(1)*y(1) x(1)*y(2) x(1)*y(3) x(2)y(1)...]
for l=1:6
    k=1;
    for i = 1:3
        for j = 1:3
            qs(l,k) = X(i,l)*Y(j,l);
            k = k+1;
        end
    end
end
Fs = null(qs);

F1 = reshape(Fs(:,1),3,3);
F2 = reshape(Fs(:,2),3,3);
F3 = reshape(Fs(:,3),3,3);

FsValues = [F1, F2, F3];
Fv1 = Fs(:,1);
Fv2 = Fs(:,2);
Fv3 = Fs(:,3);

truels = solve(Fv_true(1) - l1*Fv1(1) - l2*Fv2(1) -l3 * Fv3(1), Fv_true(2) - l1*Fv1(2) - l2*Fv2(2) -l3 * Fv3(2), Fv_true(3) - l1*Fv1(3) - l2*Fv2(3) -l3 * Fv3(3), l1,l2,l3);

l1_true = eval(truels.l1/eval(truels.l3));
l2_true = eval(truels.l2/eval(truels.l3));
F_true = F_true/eval(truels.l3);
l3_true = K_true(1,1);
% Verifying that the true values satisfy the given constraints.
res = norm(F_true - l1_true*F1 - l2_true*F2 - F3, 'fro')
paramsForTestingSparseSolver = mat2cell(reshape(FsValues,1,27), 1,ones(1,27));