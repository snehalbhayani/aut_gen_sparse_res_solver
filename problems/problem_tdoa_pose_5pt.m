function [ eqs, data0, eqs_data ] = problem_tdoa_rank2_74( data0 )

n = 5;

if nargin < 1 || isempty(data0)
    data0 = randi(50,4*n,1);
end

xx = create_vars(n);
nvect=[0 ; 0 ; 1];
U = reshape(data0(1:(3*n)),3,n);
cosangles = reshape(data0((3*n+1):(4*n)),1,n);
C = xx(1:3);
N = [xx(4:5);1];

for k = 1:n;
    eqs(k,1)=(U(:,k)-C)'*N*N'*(U(:,k)-C) - cosangles(k)^2 * (U(:,k)-C)'*(U(:,k)-C)*N'*N;
end

if nargout == 3
    xx = create_vars(n+4*n);
    oo = xx(1:n);
    data = xx((n+1):end);
    
    xx = create_vars(n);
    nvect=[0 ; 0 ; 1];
    U = reshape(data(1:(3*n)),3,n);
    cosangles = reshape(data((3*n+1):(4*n)),1,n);
    C = xx(1:3);
    N = [xx(4:5);1];
    
    for k = 1:n;
        eqs_data(k,1)=(U(:,k)-C)'*N*N'*(U(:,k)-C) - cosangles(k)^2 * (U(:,k)-C)'*(U(:,k)-C)*N'*N;
    end
end

