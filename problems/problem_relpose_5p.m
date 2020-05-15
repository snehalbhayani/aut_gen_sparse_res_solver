function [ eqs, data0, eqs_data ] = problem_relpose_5p( data0 )

if nargin < 1 || isempty(data0)
    data0 = randi(50,4*9,1);
end

xx = create_vars(3);
F0 = reshape(data0(1:9),3,3);
F1 = reshape(data0(10:18),3,3);
F2 = reshape(data0(19:27),3,3);
F3 = reshape(data0(28:36),3,3);

E = F0 + xx(1)*F1 + xx(2)*F2 +xx(3)*F3;
eqs = [2*(E*E')*E - sum(diag(E*E'))*E];
eqs = [eqs(:);det(E)];

if nargout == 3
    xx = create_vars(3+9*4);
    data = xx(4:end);

    F0 = reshape(data(1:9),3,3);
    F1 = reshape(data(10:18),3,3);
    F2 = reshape(data(19:27),3,3);
    F3 = reshape(data(28:36),3,3);

    E = F0 + xx(1)*F1 + xx(2)*F2 +xx(3)*F3;
    eqs_data = [2*(E*E')*E - sum(diag(E*E'))*E];
    eqs_data = [eqs_data(:);det(E)];
end

