function [ eqs, data0, eqs_data ] = problem_relpose_7p_fr_elimf( data0 )

if nargin < 1 || isempty(data0)
 
    data0 = randi(50,56,1);
 
end

xx = create_vars(4);
f6 = xx(1);
f7 = xx(2);
f8 = xx(3);
l = xx(4);
%w = xx(1);
f9 = 1;

h1 = data0(1:8);
h2 = data0(9:16);
h3 = data0(17:24);
h4 = data0(25:32);
h5 = data0(33:40);
h6 = data0(41:48);
h7 = data0(49:56);
vv = [l*f6 l*f7 l*f8  f6 f7 f8 l 1];
f1 = vv*h1;
f2 = vv*h2;
f3 = vv*h3;
f4 = vv*h4;
f5 = vv*h5;

F = [f1 f4 f7;f2 f5 f8;f3 f6 f9];
%K = [1 0 0;0 1 0;0 0 w];
%E = K'*F*K;

f11 = F(1,1);
f21 = F(2,1);
f31 = F(3,1);
f12 = F(1,2);
f22 = F(2,2);
f32 = F(3,2);
f13 = F(1,3);
f23 = F(2,3);
f33 = F(3,3);

%eqs = [f13*f22*f31-f12*f23*f31-f13*f21*f32+f11*f23*f32+f12*f21*f33-f11*f22*f33;...
%     f11*f13^3*f31+f13^2*f21*f23*f31+f11*f13*f23^2*f31+f21*f23^3*f31-f11*f13*f31^3-f21*f23*f31^3+f12*f13^3*f32+f13^2*f22*f23*f32+f12*f13*f23^2*f32+f22*f23^3*f32-f12*f13*f31^2*f32-f22*f23*f31^2*f32-f11*f13*f31*f32^2-f21*f23*f31*f32^2-f12*f13*f32^3-f22*f23*f32^3-f11^2*f13^2*f33-f12^2*f13^2*f33-2*f11*f13*f21*f23*f33-2*f12*f13*f22*f23*f33-f21^2*f23^2*f33-f22^2*f23^2*f33+f11^2*f31^2*f33+f21^2*f31^2*f33+2*f11*f12*f31*f32*f33+2*f21*f22*f31*f32*f33+f12^2*f32^2*f33+f22^2*f32^2*f33];

 
eqs = [f13*f22*f31-f12*f23*f31-f13*f21*f32+f11*f23*f32+f12*f21*f33-f11*f22*f33;...
      f11*f13^3*f31+f13^2*f21*f23*f31+f11*f13*f23^2*f31+f21*f23^3*f31-f11*f13*f31^3-f21*f23*f31^3+f12*f13^3*f32+f13^2*f22*f23*f32+f12*f13*f23^2*f32+f22*f23^3*f32-f12*f13*f31^2*f32-f22*f23*f31^2*f32-f11*f13*f31*f32^2-f21*f23*f31*f32^2-f12*f13*f32^3-f22*f23*f32^3-f11^2*f13^2*f33-f12^2*f13^2*f33-2*f11*f13*f21*f23*f33-2*f12*f13*f22*f23*f33-f21^2*f23^2*f33-f22^2*f23^2*f33+f11^2*f31^2*f33+f21^2*f31^2*f33+2*f11*f12*f31*f32*f33+2*f21*f22*f31*f32*f33+f12^2*f32^2*f33+f22^2*f32^2*f33];
  
 
%eqs = 2*(E*E')*E - sum(diag(E*E'))*E;
%eqs = [eqs(:);det(F)];
eqs = [eqs;l*f3-vv*h6;l^2-vv*h7];

% for iii = [1 2 4 5 9],
%     cc = coeffs(eqs(iii));
%     mm = monomials(eqs(iii));
%     mm(1,:)=mm(1,:)/2;
%     eqs(iii) = multipol(cc,mm);
% end
% for iii = [3 6 7 8],
%     cc = coeffs(eqs(iii));
%     mm = monomials(eqs(iii));
%     mm(1,:) = mm(1,:)-1;
%     mm(1,:)=mm(1,:)/2;
%     eqs(iii) = multipol(cc,mm);
% end



if nargout == 3
    xx = create_vars(4+56);
    data = xx(5:end);
    

f6 = xx(1);
f7 = xx(2);
f8 = xx(3);
l = xx(4);
%w = xx(1);
f9 = 1;

h1 = data(1:8);
h2 = data(9:16);
h3 = data(17:24);
h4 = data(25:32);
h5 = data(33:40);
h6 = data(41:48);
h7 = data(49:56);
vv = [l*f6 l*f7 l*f8  f6 f7 f8 l 1];
f1 = vv*h1;
f2 = vv*h2;
f3 = vv*h3;
f4 = vv*h4;
f5 = vv*h5;

F = [f1 f4 f7;f2 f5 f8;f3 f6 f9];
%K = [1 0 0;0 1 0;0 0 w];
%E = K'*F*K;

f11 = F(1,1);
f21 = F(2,1);
f31 = F(3,1);
f12 = F(1,2);
f22 = F(2,2);
f32 = F(3,2);
f13 = F(1,3);
f23 = F(2,3);
f33 = F(3,3);

eqs_data = [f13*f22*f31-f12*f23*f31-f13*f21*f32+f11*f23*f32+f12*f21*f33-f11*f22*f33;...
     f11*f13^3*f31+f13^2*f21*f23*f31+f11*f13*f23^2*f31+f21*f23^3*f31-f11*f13*f31^3-f21*f23*f31^3+f12*f13^3*f32+f13^2*f22*f23*f32+f12*f13*f23^2*f32+f22*f23^3*f32-f12*f13*f31^2*f32-f22*f23*f31^2*f32-f11*f13*f31*f32^2-f21*f23*f31*f32^2-f12*f13*f32^3-f22*f23*f32^3-f11^2*f13^2*f33-f12^2*f13^2*f33-2*f11*f13*f21*f23*f33-2*f12*f13*f22*f23*f33-f21^2*f23^2*f33-f22^2*f23^2*f33+f11^2*f31^2*f33+f21^2*f31^2*f33+2*f11*f12*f31*f32*f33+2*f21*f22*f31*f32*f33+f12^2*f32^2*f33+f22^2*f32^2*f33];

%eqs = 2*(E*E')*E - sum(diag(E*E'))*E;
%eqs = [eqs(:);det(F)];
eqs_data = [eqs_data;l*f3-vv*h6;l^2-vv*h7];

% for iii = [1 2 4 5 9],
%     cc = coeffs(eqs(iii));
%     mm = monomials(eqs(iii));
%     mm(1,:)=mm(1,:)/2;
%     eqs(iii) = multipol(cc,mm);
% end
% for iii = [3 6 7 8],
%     cc = coeffs(eqs(iii));
%     mm = monomials(eqs(iii));
%     mm(1,:) = mm(1,:)-1;
%     mm(1,:)=mm(1,:)/2;
%     eqs(iii) = multipol(cc,mm);
% end


    totalsyms = 4 + 56;
    noofvars = 4;
    ipparams = strcat('a',num2str(1));
    for i = 2:totalsyms
        if i > noofvars
            ipparams = strjoin({ipparams, ',c', num2str(i-noofvars)}, '');
        else
            ipparams = strjoin({ipparams, ',a', num2str(i)}, '');
        end
    end
    
    fileID = fopen('Eqs_problem_relpose_7p_fr_elimf.m','w');
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

