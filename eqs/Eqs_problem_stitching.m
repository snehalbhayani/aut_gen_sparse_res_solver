function eqs = retrieve_eqs(a1,a2,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,c31,c32) 
eqs(1) = a1^6*a2^3*c1 + a1^5*a2^3*c2 + a1^4*a2^3*c3 + a1^4*a2^2*c4 + a1^3*a2^3*c5 + a1^3*a2^2*c6 + a1^2*a2^3*c7 + a1^2*a2^2*c8 + a1^2*a2*c9 + a1*a2^3*c10 + a1*a2^2*c11 + a1*a2*c12 + a2^3*c13 + a2^2*c14 + a2*c15 + c16;
eqs(2) = a1^6*a2^3*c17 + a1^5*a2^3*c18 + a1^4*a2^3*c19 + a1^4*a2^2*c20 + a1^3*a2^3*c21 + a1^3*a2^2*c22 + a1^2*a2^3*c23 + a1^2*a2^2*c24 + a1^2*a2*c25 + a1*a2^3*c26 + a1*a2^2*c27 + a1*a2*c28 + a2^3*c29 + a2^2*c30 + a2*c31 + c32;

end