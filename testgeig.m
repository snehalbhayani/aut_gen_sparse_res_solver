%%
clear all;
clc;
% Size of the matrices to be tested for generalized eigenvalues.
nge = 40;
% Size of the matrices to be tested for normal eigenvalues.
ne = 40;
% Number of matrices to test for or number of iterations of the sam matrix
% to test for.
iters = 1000;
fprintf("Testing for random %d iterations for matrices of dimensions %d* %d \n", iters, nge, nge);
%% generate 1000 random matrices
% Create a random vector of size nge * nge elements and we then reshape it
% to dimensions, nge X nge. 
A = randn(nge^2, iters);
% A = repmat(randn(nge^2,1),1,iters);
A = transpose(mat2cell(reshape(A, nge, nge*iters), nge, nge*ones(1,iters)));

% Create a random vector of size nge * nge elements and we then reshape it
% to dimensions, nge X nge. 
B = randn(nge^2, iters);
% B = repmat(randn(nge^2,1),1,iters);
B = transpose(mat2cell(reshape(B, nge, nge*iters), nge, nge*ones(1,iters)));

% Create a random vector of size ne * ne elements and we then reshape it
% to dimensions, ne X ne. 
C = randn(ne^2, iters);
% C = repmat(randn(ne^2,1),1,iters);
C = transpose(mat2cell(reshape(C, ne, ne*iters), ne, ne*ones(1,iters)));

ta = A{1};
tb = B{1};


for i = 1:size(A,1)
    tic;
    for j = 1:i
        [V,D] = eig(A{j}, B{j});
    end
    tg(i) = toc/i;
end

tc = C{1};
for i = 1:size(C,1)
    tic;
    for j = 1:i
        [V,D] = eig(C{j});
    end
    te(i) = toc/i;

end

%% Plots
rang = 100:iters-100;
figure(1);
clf;
plot(rang, tg(rang), rang, te(rang));
title("Plot of time taken for eig");
xlabel('iter');
ylabel('time taken in secs');
legend({'Gen. eig.','Eig.'},'Location','southwest');

