function [ eqs, data0, eqs_data ] = problem_relpose_7p_fr_1s_elim( data0 )

if nargin < 1 || isempty(data0)
    data0 = randi(50,5*12,1);
end

xx = create_vars(4);

F0 = reshape(data0(1:12),3,4);
F1 = reshape(data0(13:24),3,4);
F2 = reshape(data0(25:36),3,4);
F3 = reshape(data0(37:48),3,4);
F4 = reshape(data0(49:60),3,4);

F = F0 + xx(1)*F1 + xx(2)*F2 + xx(3)*F3 + xx(4)*F4;

f11 = F(1,1);
f21 = F(2,1);
f31 = F(3,1);
f12 = F(1,2);
f22 = F(2,2);
f32 = F(3,2);
f13 = F(1,3);
f23 = F(2,3);
f33 = F(3,3);
y13 = F(1,4);
y23 = F(2,4);
y33 = F(3,4);

eqs = [f33*y23-f23*y33;...
      f33*y13-f13*y33;...
      f23*y13-f13*y23;...
      f22*f31*y13-f21*f32*y13-f12*f31*y23+f11*f32*y23+f12*f21*y33-f11*f22*y33;...
      f13*f22*f31-f12*f23*f31-f13*f21*f32+f11*f23*f32+f12*f21*f33-f11*f22*f33;...
      f11*f31*y13*y23+f12*f32*y13*y23+f21*f31*y23^2+f22*f32*y23^2-f11*f21*y13*y33-f12*f22*y13*y33-f21^2*y23*y33-f22^2*y23*y33+f31^2*y23*y33+f32^2*y23*y33-f21*f31*y33^2-f22*f32*y33^2;...
      f11*f13*f31*y23+f21*f23*f31*y23+f12*f13*f32*y23+f22*f23*f32*y23-f11*f13*f21*y33-f12*f13*f22*y33-f21^2*f23*y33-f22^2*f23*y33+f23*f31^2*y33+f23*f32^2*y33-f21*f31*f33*y33-f22*f32*f33*y33;...
      f11*f31*y13^2+f12*f32*y13^2+f21*f31*y13*y23+f22*f32*y13*y23-f11^2*y13*y33-f12^2*y13*y33+f31^2*y13*y33+f32^2*y13*y33-f11*f21*y23*y33-f12*f22*y23*y33-f11*f31*y33^2-f12*f32*y33^2;...
      f11*f21*y13^2+f12*f22*y13^2-f11^2*y13*y23-f12^2*y13*y23+f21^2*y13*y23+f22^2*y13*y23-f11*f21*y23^2-f12*f22*y23^2+f21*f31*y13*y33+f22*f32*y13*y33-f11*f31*y23*y33-f12*f32*y23*y33;...
      f11*f13*f31*y13+f12*f13*f32*y13+f13*f21*f31*y23+f13*f22*f32*y23-f11^2*f13*y33-f12^2*f13*y33-f11*f21*f23*y33-f12*f22*f23*y33+f13*f31^2*y33+f13*f32^2*y33-f11*f31*f33*y33-f12*f32*f33*y33;...
      f11*f13*f21*y13+f12*f13*f22*y13-f11^2*f13*y23-f12^2*f13*y23+f13*f21^2*y23+f13*f22^2*y23-f11*f21*f23*y23-f12*f22*f23*y23+f13*f21*f31*y33-f11*f23*f31*y33+f13*f22*f32*y33-f12*f23*f32*y33;...
      f11*f13*f23*f31+f21*f23^2*f31+f12*f13*f23*f32+f22*f23^2*f32-f11*f13*f21*f33-f12*f13*f22*f33-f21^2*f23*f33-f22^2*f23*f33+f23*f31^2*f33+f23*f32^2*f33-f21*f31*f33^2-f22*f32*f33^2;...
      f11*f13^2*f31+f13*f21*f23*f31+f12*f13^2*f32+f13*f22*f23*f32-f11^2*f13*f33-f12^2*f13*f33-f11*f21*f23*f33-f12*f22*f23*f33+f13*f31^2*f33+f13*f32^2*f33-f11*f31*f33^2-f12*f32*f33^2;...
      f11*f13^2*f21+f12*f13^2*f22-f11^2*f13*f23-f12^2*f13*f23+f13*f21^2*f23+f13*f22^2*f23-f11*f21*f23^2-f12*f22*f23^2+f13*f21*f31*f33-f11*f23*f31*f33+f13*f22*f32*f33-f12*f23*f32*f33];
  

if nargout == 3
    xx = create_vars(4+60);
    data = xx(5:end);
    
 
F0 = reshape(data(1:12),3,4);
F1 = reshape(data(13:24),3,4);
F2 = reshape(data(25:36),3,4);
F3 = reshape(data(37:48),3,4);
F4 = reshape(data(49:60),3,4);

F = F0 + xx(1)*F1 + xx(2)*F2 + xx(3)*F3 + xx(4)*F4;

f11 = F(1,1);
f21 = F(2,1);
f31 = F(3,1);
f12 = F(1,2);
f22 = F(2,2);
f32 = F(3,2);
f13 = F(1,3);
f23 = F(2,3);
f33 = F(3,3);
y13 = F(1,4);
y23 = F(2,4);
y33 = F(3,4);

eqs_data = [f33*y23-f23*y33;...
      f33*y13-f13*y33;...
      f23*y13-f13*y23;...
      f22*f31*y13-f21*f32*y13-f12*f31*y23+f11*f32*y23+f12*f21*y33-f11*f22*y33;...
      f13*f22*f31-f12*f23*f31-f13*f21*f32+f11*f23*f32+f12*f21*f33-f11*f22*f33;...
      f11*f31*y13*y23+f12*f32*y13*y23+f21*f31*y23^2+f22*f32*y23^2-f11*f21*y13*y33-f12*f22*y13*y33-f21^2*y23*y33-f22^2*y23*y33+f31^2*y23*y33+f32^2*y23*y33-f21*f31*y33^2-f22*f32*y33^2;...
      f11*f13*f31*y23+f21*f23*f31*y23+f12*f13*f32*y23+f22*f23*f32*y23-f11*f13*f21*y33-f12*f13*f22*y33-f21^2*f23*y33-f22^2*f23*y33+f23*f31^2*y33+f23*f32^2*y33-f21*f31*f33*y33-f22*f32*f33*y33;...
      f11*f31*y13^2+f12*f32*y13^2+f21*f31*y13*y23+f22*f32*y13*y23-f11^2*y13*y33-f12^2*y13*y33+f31^2*y13*y33+f32^2*y13*y33-f11*f21*y23*y33-f12*f22*y23*y33-f11*f31*y33^2-f12*f32*y33^2;...
      f11*f21*y13^2+f12*f22*y13^2-f11^2*y13*y23-f12^2*y13*y23+f21^2*y13*y23+f22^2*y13*y23-f11*f21*y23^2-f12*f22*y23^2+f21*f31*y13*y33+f22*f32*y13*y33-f11*f31*y23*y33-f12*f32*y23*y33;...
      f11*f13*f31*y13+f12*f13*f32*y13+f13*f21*f31*y23+f13*f22*f32*y23-f11^2*f13*y33-f12^2*f13*y33-f11*f21*f23*y33-f12*f22*f23*y33+f13*f31^2*y33+f13*f32^2*y33-f11*f31*f33*y33-f12*f32*f33*y33;...
      f11*f13*f21*y13+f12*f13*f22*y13-f11^2*f13*y23-f12^2*f13*y23+f13*f21^2*y23+f13*f22^2*y23-f11*f21*f23*y23-f12*f22*f23*y23+f13*f21*f31*y33-f11*f23*f31*y33+f13*f22*f32*y33-f12*f23*f32*y33;...
      f11*f13*f23*f31+f21*f23^2*f31+f12*f13*f23*f32+f22*f23^2*f32-f11*f13*f21*f33-f12*f13*f22*f33-f21^2*f23*f33-f22^2*f23*f33+f23*f31^2*f33+f23*f32^2*f33-f21*f31*f33^2-f22*f32*f33^2;...
      f11*f13^2*f31+f13*f21*f23*f31+f12*f13^2*f32+f13*f22*f23*f32-f11^2*f13*f33-f12^2*f13*f33-f11*f21*f23*f33-f12*f22*f23*f33+f13*f31^2*f33+f13*f32^2*f33-f11*f31*f33^2-f12*f32*f33^2;...
      f11*f13^2*f21+f12*f13^2*f22-f11^2*f13*f23-f12^2*f13*f23+f13*f21^2*f23+f13*f22^2*f23-f11*f21*f23^2-f12*f22*f23^2+f13*f21*f31*f33-f11*f23*f31*f33+f13*f22*f32*f33-f12*f23*f32*f33];
  

    totalsyms = 4 + 60;
    noofvars = 4;
    ipparams = strcat('a',num2str(1));
    for i = 2:totalsyms
        if i > noofvars
            ipparams = strjoin({ipparams, ',c', num2str(i-noofvars)}, '');
        else
            ipparams = strjoin({ipparams, ',a', num2str(i)}, '');
        end
    end
    
    fileID = fopen('Eqs_problem_relpose_7p_fr_1s_elim.m','w');
    fprintf(fileID, '%s', strjoin({'function eqs = retrieve_eqs(', ipparams, ') '},''));
    fprintf(fileID, '\n');
    for i = 1:size(eqs,1)
        fprintf(fileID, '%s', strjoin({'eqs(', num2str(i),') = ',char(eqs_data(i,1), true, [], true, noofvars), ';'},''));
        fprintf(fileID, '\n');
    end
    fprintf(fileID, '\n');
    fprintf(fileID, '%s', 'end');
    fclose(fileID);

end

