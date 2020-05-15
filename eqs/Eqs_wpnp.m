function eqs = retrieve_eqs(a1,a2,a3,a4,c1,c2,c3,c4,c5,c6,c7,c8,c9)
eqs(1) = 4*a1^3*c1^2 + 4*a1^3*c2^2 + 4*a1*a2^2*c1^2 - 4*a1*a2^2*c2^2 + 8*a1*a2^2*c3^2 - 4*a1*a3^2*c1^2 + 4*a1*a3^2*c2^2 + 8*a1*a3^2*c3^2 + 4*a1*a4^2*c1^2 + 4*a1*a4^2*c2^2 - 4*a1*c1*c4 - 4*a1*c2*c7 + 8*a2*a3*a4*c1^2 - 8*a2*a3*a4*c2^2 + 4*a2*c3*c9 - 4*a3*c3*c8 - 4*a4*c1*c5 + 4*a4*c2*c6;eqs(2) = 4*a1^2*a2*c1^2 - 4*a1^2*a2*c2^2 + 8*a1^2*a2*c3^2 + 8*a1*a3*a4*c1^2 - 8*a1*a3*a4*c2^2 + 4*a1*c3*c9 + 4*a2^3*c1^2 + 4*a2^3*c2^2 + 4*a2*a3^2*c1^2 + 4*a2*a3^2*c2^2 - 4*a2*a4^2*c1^2 + 4*a2*a4^2*c2^2 + 8*a2*a4^2*c3^2 - 4*a2*c1*c4 + 4*a2*c2*c7 - 4*a3*c1*c5 - 4*a3*c2*c6 - 4*a4*c3*c8;eqs(3) = -4*a1^2*a3*c1^2 + 4*a1^2*a3*c2^2 + 8*a1^2*a3*c3^2 + 8*a1*a2*a4*c1^2 - 8*a1*a2*a4*c2^2 - 4*a1*c3*c8 + 4*a2^2*a3*c1^2 + 4*a2^2*a3*c2^2 - 4*a2*c1*c5 - 4*a2*c2*c6 + 4*a3^3*c1^2 + 4*a3^3*c2^2 + 4*a3*a4^2*c1^2 - 4*a3*a4^2*c2^2 + 8*a3*a4^2*c3^2 + 4*a3*c1*c4 - 4*a3*c2*c7 - 4*a4*c3*c9;eqs(4) = 4*a1^2*a4*c1^2 + 4*a1^2*a4*c2^2 + 8*a1*a2*a3*c1^2 - 8*a1*a2*a3*c2^2 - 4*a1*c1*c5 + 4*a1*c2*c6 - 4*a2^2*a4*c1^2 + 4*a2^2*a4*c2^2 + 8*a2^2*a4*c3^2 - 4*a2*c3*c8 + 4*a3^2*a4*c1^2 - 4*a3^2*a4*c2^2 + 8*a3^2*a4*c3^2 - 4*a3*c3*c9 + 4*a4^3*c1^2 + 4*a4^3*c2^2 + 4*a4*c1*c4 + 4*a4*c2*c7;
end