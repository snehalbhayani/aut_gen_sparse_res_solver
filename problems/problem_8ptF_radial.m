function [ eqs, data0, eqs_data ] = problem_8ptF_radial( data0 )

if nargin < 1 || isempty(data0)
    data0 = randi(50,50,1);
end

d1 = data0(1:7);
d2 = data0(7+(1:7));
d3 = data0(14+(1:7));
d4 = data0(21+(1:7));
d5 = data0(28+(1:7));
d6 = data0(35+(1:7));

xx = create_vars(3);

lam = xx(3);

f11 = d1(1)*xx(1)*xx(3)+d1(2)*xx(2)*xx(3)+d1(3)*xx(3)^2+d1(4)*xx(1)+d1(5)*xx(2)+d1(6)*xx(3)+d1(7);
f21 = d2(1)*xx(1)*xx(3)+d2(2)*xx(2)*xx(3)+d2(3)*xx(3)^2+d2(4)*xx(1)+d2(5)*xx(2)+d2(6)*xx(3)+d2(7);
f31 = d3(1)*xx(1)*xx(3)+d3(2)*xx(2)*xx(3)+d3(3)*xx(3)^2+d3(4)*xx(1)+d3(5)*xx(2)+d3(6)*xx(3)+d3(7);
f12 = d4(1)*xx(1)*xx(3)+d4(2)*xx(2)*xx(3)+d4(3)*xx(3)^2+d4(4)*xx(1)+d4(5)*xx(2)+d4(6)*xx(3)+d4(7);
f22 = d5(1)*xx(1)*xx(3)+d5(2)*xx(2)*xx(3)+d5(3)*xx(3)^2+d5(4)*xx(1)+d5(5)*xx(2)+d5(6)*xx(3)+d5(7);
f32 = d6(1)*xx(1)*xx(3)+d6(2)*xx(2)*xx(3)+d6(3)*xx(3)^2+d6(4)*xx(1)+d6(5)*xx(2)+d6(6)*xx(3)+d6(7);

f13 = xx(1);
f23 = xx(2);
f33 = 1;


x1 = reshape(data0(43:46),2,2);
x2 = reshape(data0(47:50),2,2);

x1 = [x1;1+(x1(1,:).^2+x1(2,:).^2)*lam];
x2 = [x2;1+(x2(1,:).^2+x2(2,:).^2)*lam];
eqs(3,1)=det([f11 f12 f13;f21 f22 f23;f31 f32 f33]);
fv = [f11;f21;f31;f12;f22;f32;f13;f23;f33];
for iii = 1:2,
    tmp = x1(:,iii)*x2(:,iii)';
    eqs(iii)=tmp(:)'*fv;
end



if nargout == 3
    xx = create_vars(3+50);
    data = xx(4:end);
    
  
d1 = data(1:7);
d2 = data(7+(1:7));
d3 = data(14+(1:7));
d4 = data(21+(1:7));
d5 = data(28+(1:7));
d6 = data(35+(1:7));


lam = xx(3);

f11 = d1(1)*xx(1)*xx(3)+d1(2)*xx(2)*xx(3)+d1(3)*xx(3)^2+d1(4)*xx(1)+d1(5)*xx(2)+d1(6)*xx(3)+d1(7);
f21 = d2(1)*xx(1)*xx(3)+d2(2)*xx(2)*xx(3)+d2(3)*xx(3)^2+d2(4)*xx(1)+d2(5)*xx(2)+d2(6)*xx(3)+d2(7);
f31 = d3(1)*xx(1)*xx(3)+d3(2)*xx(2)*xx(3)+d3(3)*xx(3)^2+d3(4)*xx(1)+d3(5)*xx(2)+d3(6)*xx(3)+d3(7);
f12 = d4(1)*xx(1)*xx(3)+d4(2)*xx(2)*xx(3)+d4(3)*xx(3)^2+d4(4)*xx(1)+d4(5)*xx(2)+d4(6)*xx(3)+d4(7);
f22 = d5(1)*xx(1)*xx(3)+d5(2)*xx(2)*xx(3)+d5(3)*xx(3)^2+d5(4)*xx(1)+d5(5)*xx(2)+d5(6)*xx(3)+d5(7);
f32 = d6(1)*xx(1)*xx(3)+d6(2)*xx(2)*xx(3)+d6(3)*xx(3)^2+d6(4)*xx(1)+d6(5)*xx(2)+d6(6)*xx(3)+d6(7);

f13 = xx(1);
f23 = xx(2);
f33 = 1;


x1 = reshape(data(43:46),2,2);
x2 = reshape(data(47:50),2,2);

x1 = [x1;1+(x1(1,:).^2+x1(2,:).^2)*lam];
x2 = [x2;1+(x2(1,:).^2+x2(2,:).^2)*lam];
eqs_data(3,1)=det([f11 f12 f13;f21 f22 f23;f31 f32 f33]);
fv = [f11;f21;f31;f12;f22;f32;f13;f23;f33];
for iii = 1:2,
    tmp = x1(:,iii)*x2(:,iii)';
    eqs_data(iii)=tmp(:)'*fv;
end




end

