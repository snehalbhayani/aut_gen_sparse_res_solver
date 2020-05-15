%% instructions for benchmarking of the time performance for the a solver using Sparse Resultant based solver.
% As against the benchmarking of Groebner basis algorithm, this
% benchmarking code is generic in nature and hence can be used to benchmark
% a solver for any problem.

% parameters : 1. As -- A vertical stack of the A matrices, of the form [A1;A2;...]
% parameters : 2. Bs -- A vertical stack of the B matrices of the form [B1;
% B2; ... ]; Also if the number of iterations is iters, then we have
% size(As,1) = size(As,2) * iters;

% The generalized eigen problem being solved is of the form, A1*v = lam * B1 * v;

% The function prints, the total time taken for eig. solutions of iters
% runs and also the mean time. 


function[] = benchmark_solver(As,Bs)
%%

i1 = 1; 
i2 = 2;
j1 = 2;
j2 = 1;
thres = 1e-3;

allAs = mat2cell(As,size(As,2)*ones(1,size(As,1)/size(As,2)), [size(As,2)]);
allBs = mat2cell(Bs,size(Bs,2)*ones(1,size(Bs,1)/size(Bs,2)), [size(Bs,2)]);
tic;
for i=1:size(allAs,1)
    for j = 1:size(allAs,2)
        [V,D] = eig(allAs{i,1},allBs{i,1});
%         V = V ./V(1,:);
        V = V(:,find(abs(V(i2,:)./V(i1,:) - V(j2,:)./V(j1,:)) <= thres));
    end
end
elasped_time = toc;
disp("Time taken for all iters is ::: "+ elasped_time);
disp("Time taken per iter is ::: "+ elasped_time/size(allAs,1));


end