function [u, X, P] = find_v(n, r, C, c13_t, c14_t, thresh, d)
u = randn(3,1); u = u/norm(u);
ntu = (transpose(n)*u);
P = C - u * (transpose(n)*C + d)/ntu;

syms vx vy vz;
tht1 = abs(acos(dot(n, u)/(norm(n))) - pi);
tht2 = asin(sin(tht1)*r);
v_sym = [vx;vy;vz];
eq1 = (transpose(n) * v_sym)^2 - (norm(n)*cos(tht2))^2;
    eq2 = transpose(cross(u,n)) * v_sym;
    eq3 = vx^2 + vy^2 + vz^2 - 1;
    [x,y,z] = solve([eq1;eq2;eq3], [vx;vy;vz]);
    for i=1:length(x)
        v = eval([x(i);y(i);z(i)]);
        X = P - rand() * v;
        if (transpose(n)*C+d)*(transpose(n)*X+d) < 0 && norm(r * cross(u,n) - cross(v,n)) < thresh
            break;
        else
            X = [];
            v = [];
        end
    end
end