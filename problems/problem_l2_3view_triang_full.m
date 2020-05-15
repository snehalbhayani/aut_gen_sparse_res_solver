function [ eqs, data0, eqs_data ] = problem_l2_3view_triang_full( data0 )

if nargin < 1 || isempty(data0)
    
%     rr = [3 -3 4 -4 5 -5 12 -12 7 -7 24 -24;4 4 3 3 12 12 5 5 24 24 7 7;5 5 5 5 13 13 13 13 25 25 25 25  ]; 
%     ii = randi(12,1,6);
%     rc = rr(1,ii);
%     rs = rr(2,ii);
%     rw = rr(3,ii);
%     R1 = [rc(1) rs(1) 0;-rs(1) rc(1) 0;0 0 rw(1)]*[rc(2) 0 rs(2);0 rw(2) 0;-rs(2) 0 rc(2)]*[rw(3) 0 0;0 rc(3) rs(3);0 -rs(3) rc(3)];
%     R2 = [rc(4) rs(4) 0;-rs(4) rc(4) 0;0 0 rw(4)]*[rc(5) 0 rs(5);0 rw(5) 0;-rs(5) 0 rc(5)]*[rw(6) 0 0;0 rc(6) rs(6);0 -rs(6) rc(6)];
%     R3 = R2*R1;
%     t1 = randi(50,3,1);
%     t2 = randi(50,3,1);
%     t3 = R2*t1+t2;
%    
%     E12 = R1*[0 -t1(3) t1(2);t1(3) 0 -t1(1);-t1(2) t1(1) 0];
%     E23 = R2*[0 -t2(3) t2(2);t2(3) 0 -t2(1);-t2(2) t2(1) 0];
%     E13 = R3*[0 -t3(3) t3(2);t3(3) 0 -t3(1);-t3(2) t3(1) 0];
      E12 = randi(20,3,3);
      E23 = randi(20,3,3);
      E13 = randi(20,3,3);
        
     data0 = [E12(:);E23(:);E13(:);randi(50,6,1)];
end

E12 = reshape(data0(1:9),3,3);
E23 = reshape(data0(10:18),3,3);
E13 = reshape(data0(19:27),3,3);
u = [reshape(data0(28:33),2,3);ones(1,3)];

xx = create_vars(9);
uh = [reshape(xx(1:6),2,3);ones(1,3)];
l1 = xx(7);
l2 = xx(8);
l3 = xx(9);
% S = [1 0 0;0 1 0];
% 
% eqs = [uh(:,1)'*E12*uh(:,2);uh(:,2)'*E23*uh(:,3)];
% eqs = [eqs;2*S*(uh(:,1)-u(:,1))+l1*S*E12*uh(:,2)];
% eqs = [eqs;2*S*(uh(:,2)-u(:,2))+l1*S*E12'*uh(:,1)+l2*S*E23*uh(:,3)];
% eqs = [eqs;2*S*(uh(:,3)-u(:,3))+l2*S*E23'*uh(:,2)];
f = sum((uh(:)-u(:)).^2)+l1*uh(:,1)'*E12*uh(:,2)+l2*uh(:,2)'*E23*uh(:,3)+l3*uh(:,1)'*E13*uh(:,3);
eqs = diff(f);
eqs = eqs(:);

if nargout == 3
    xx = create_vars(9+27+6);
    data = xx(10:end);
    
E12 = reshape(data(1:9),3,3);
E23 = reshape(data(10:18),3,3);
E13 = reshape(data(19:27),3,3);
u = [reshape(data(28:33),2,3);ones(1,3)];
uh = [reshape(xx(1:6),2,3);ones(1,3)];
l1 = xx(7);
l2 = xx(8);
l3 = xx(9);
% S = [1 0 0;0 1 0];
% 
% eqs_data = [uh(:,1)'*E12*uh(:,2);uh(:,2)'*E23*uh(:,3)];
% eqs_data = [eqs_data;2*S*(uh(:,1)-u(:,1))+l1*S*E12*uh(:,2)];
% eqs_data = [eqs_data;2*S*(uh(:,2)-u(:,2))+l1*S*E12'*uh(:,1)+l2*S*E23*uh(:,3)];
% eqs_data = [eqs_data;2*S*(uh(:,3)-u(:,3))+l2*S*E23'*uh(:,2)];
f = sum((uh(:)-u(:)).^2)+l1*uh(:,1)'*E12*uh(:,2)+l2*uh(:,2)'*E23*uh(:,3)+l3*uh(:,1)'*E13*uh(:,3);
eqs_data = diff(f);
eqs_data = eqs_data(:);



end

