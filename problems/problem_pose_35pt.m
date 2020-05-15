function [ eqs, data0, eqs_data ] = problem_pose_35pt( data0 )

if nargin < 1 || isempty(data0)
    
    data0 = randi(50,4*3+3.5*2,1);
 
end

X = reshape(data0(1:12),3,4);
x = data0(13:16);
y = data0(17:19);

xx = create_vars(2);
qx = xx(1);
qy = xx(2);

R = [1+qx^2-qy^2 2*qx*qy 2*qy;...
    2*qx*qy qy^2-qx^2+1 -2*qx;...
    -2*qy 2*qx 1-qx^2-qy^2];


R1 = R(1,:);
R2 = R(2,:);
R3 = R(3,:);
F(4,3)=multipol;

iddes = [1 2 3;1 3 2;2 4 3;3 4 2];

for lll = 1:4,
    i = iddes(lll,1);
    j = iddes(lll,2);
    k = iddes(lll,3);
    

    F(lll,1)=(y(i)-y(k))*R1*(X(:,j)-X(:,i))-(x(i)-x(j))*R2*(X(:,k)-X(:,i));
    F(lll,2)=-(y(i)-y(k))*R2*(X(:,j)-X(:,i))-(x(i)-x(j))*R1*(X(:,k)-X(:,i));
    F(lll,3)=-(y(i)-y(k))*x(j)*R3*(X(:,j)-X(:,i))+(x(i)-x(j))*y(k)*R3*(X(:,k)-X(:,i));
end
eqs = [];
for iii = 1:4,
    eqs = [eqs;det(F(setdiff(1:4,iii),:))];
end



if nargout == 3
    xx = create_vars(2+19);
    data = xx(3:end);
    
X = reshape(data(1:12),3,4);
x = data(13:16);
y = data(17:19);


qx = xx(1);
qy = xx(2);

R = [1+qx^2-qy^2 2*qx*qy 2*qy;...
    2*qx*qy qy^2-qx^2+1 -2*qx;...
    -2*qy 2*qx 1-qx^2-qy^2];


R1 = R(1,:);
R2 = R(2,:);
R3 = R(3,:);
F(4,3)=multipol;

iddes = [1 2 3;1 3 2;2 4 3;3 4 2];

for lll = 1:4,
    i = iddes(lll,1);
    j = iddes(lll,2);
    k = iddes(lll,3);
    

    F(lll,1)=(y(i)-y(k))*R1*(X(:,j)-X(:,i))-(x(i)-x(j))*R2*(X(:,k)-X(:,i));
    F(lll,2)=-(y(i)-y(k))*R2*(X(:,j)-X(:,i))-(x(i)-x(j))*R1*(X(:,k)-X(:,i));
    F(lll,3)=-(y(i)-y(k))*x(j)*R3*(X(:,j)-X(:,i))+(x(i)-x(j))*y(k)*R3*(X(:,k)-X(:,i));
end
eqs_data = [];
for iii = 1:4,
    disp(iii)
    eqs_data = [eqs_data;det(F(setdiff(1:4,iii),:))];
end



end

