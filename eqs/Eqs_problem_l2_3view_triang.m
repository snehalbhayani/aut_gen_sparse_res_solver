function eqs = retrieve_eqs(a1,a2,a3,a4,a5,a6,a7,a8,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24) 
eqs(1) = 2*a1 + a3*a7*c1 + a4*a7*c4 + a7*c7 - 2*c19;
eqs(2) = 2*a2 + a3*a7*c2 + a4*a7*c5 + a7*c8 - 2*c20;
eqs(3) = a1*a7*c1 + a2*a7*c2 + 2*a3 + a5*a8*c10 + a6*a8*c13 + a7*c3 + a8*c16 - 2*c21;
eqs(4) = a1*a7*c4 + a2*a7*c5 + 2*a4 + a5*a8*c11 + a6*a8*c14 + a7*c6 + a8*c17 - 2*c22;
eqs(5) = a3*a8*c10 + a4*a8*c11 + 2*a5 + a8*c12 - 2*c23;
eqs(6) = a3*a8*c13 + a4*a8*c14 + 2*a6 + a8*c15 - 2*c24;
eqs(7) = a1*a3*c1 + a1*a4*c4 + a1*c7 + a2*a3*c2 + a2*a4*c5 + a2*c8 + a3*c3 + a4*c6 + c9;
eqs(8) = a3*a5*c10 + a3*a6*c13 + a3*c16 + a4*a5*c11 + a4*a6*c14 + a4*c17 + a5*c12 + a6*c15 + c18;

end