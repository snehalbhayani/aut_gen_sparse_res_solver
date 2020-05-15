function [eqs,data0,eqs_data] = problem_example2_hv(data0)
% Same as example 2, but we have eliminated z using hidden variable.

if nargin < 1 || isempty(data0)
    % no input, generate a random integer instance
    data0 = randi(100,3*4*3,1);
end

% Setup equation system
Q1 = reshape(data0(1:12),3,4);
Q2 = reshape(data0(13:24),3,4);
Q3 = reshape(data0(25:36),3,4);

xx = create_vars(3);
eqs_orig(1) = [xx(1:2);1]' * Q1 * [xx;1];
eqs_orig(2) = [xx(1:2);1]' * Q2 * [xx;1];
eqs_orig(3) = [xx(1:2);1]' * Q3 * [xx;1];

% Hide variable z.
% Rewrite equations as A*[z;1] = 0 where A depends on x and y
A = collect_terms(eqs_orig, [xx(3); 1]);
% and then require all sub-determinants of A to vanish
ind = nchoosek(1:size(A,1),size(A,2));
eqs = [];
for k = 1:size(ind,1)    
    eqs = [eqs; det(A(ind(k,:),:))];
    % If you have trouble with integer overflow you can use
    % eqs = [eqs; zp_det(A(ind(k,:),:),30097)];
end
eqs = remove_vars(eqs,3); % remove z from our equations

% Setup equation with data as additional unknowns
if nargout == 3
    xx0 = create_vars(3+3*4*3);
    xx = xx0(1:3);
    data = xx0(4:end);

    % Setup equation system
    Q1 = reshape(data(1:12),3,4);
    Q2 = reshape(data(13:24),3,4);
    Q3 = reshape(data(25:36),3,4);

    xx = create_vars(3);
    eqs_orig(1) = [xx(1:2);1]' * Q1 * [xx;1];
    eqs_orig(2) = [xx(1:2);1]' * Q2 * [xx;1];
    eqs_orig(3) = [xx(1:2);1]' * Q3 * [xx;1];

    % Hide variable z
    A = collect_terms(eqs_orig, [xx(3); 1]);
    ind = nchoosek(1:size(A,1),size(A,2))
    eqs_data = [];
    for k = 1:size(ind,1)
        eqs_data = [eqs_data; det(A(ind(k,:),:))];
    end
    eqs_data = remove_vars(eqs_data,3); % remove z from our equations

end