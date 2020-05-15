function [ eqs, data0, eqs_data ] = problem_satellite_triang( data0 )

if nargin < 1 || isempty(data0)
    data0 = randi(100,20*3,1);
end


%xx = create_vars(3);
mm = [0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 2 2 2 3;...
    0 0 0 0 1 1 1 2 2 3 0 0 0 1 1 2 0 0 1 0 ;...
    0 1 2 3 0 1 2 0 1 0 0 1 2 0 1 0 0 1 0 0 ];



c1 = data0(1:20);
c2 = data0(21:40);
c3 = data0(41:60);

eqs = [multipol(c1',mm);multipol(c2',mm);multipol(c3',mm)];

if nargout == 3
    xx = create_vars(3+60);
    data = xx(4:end);
    c1 = data(1:20);
    c2 = data(21:40);
    c3 = data(41:60);
    eqs_data = [sum(c1'.*xx(1).^mm(1,:).*xx(2).^mm(2,:).*xx(3).^mm(3,:));...
        sum(c2'.*xx(1).^mm(1,:).*xx(2).^mm(2,:).*xx(3).^mm(3,:));...
        sum(c3'.*xx(1).^mm(1,:).*xx(2).^mm(2,:).*xx(3).^mm(3,:))];
end
