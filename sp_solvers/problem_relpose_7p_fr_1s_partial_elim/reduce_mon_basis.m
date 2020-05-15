function[indicesToRemove, numOfRedCoeffMat, zeigvalindx, solverTemplate, extendedbasis, nullspacesize1, depdXcols1, nzrows1, depdCcols1, indepdCcols1, mdepdXcols1, mindepdXcols1, indepdCrows1, rowstorem1, nullspacesize2, depdXcols2, nzrows2, depdCcols2, indepdCcols2, mdepdXcols2, mindepdXcols2, indepdCrows2, rowstorem2, nullspacesize3, depdXcols3, nzrows3, depdCcols3, indepdCcols3, mdepdXcols3, mindepdXcols3, indepdCrows3, rowstorem3] = reducemonbasis(data, noofdatacoeff, foldername)
hiddenvarnumber = 4;
for i = 1:length(data)
eval(strjoin({'c',num2str(i),' = ', 'data(',num2str(i),');'},''));
end
t1 = c2;
t2 = c1;
t3 = c3;
t4 = c4;
t5 = c5;
t6 = -c10*c30+c15*c25;
t7 = -c10*c28+c10*c35+c13*c25+c15*c23-c20*c25-c30*c8;
t8 = -c10*c29+c14*c25+c15*c24-c30*c9;
t9 = -c10*c27+c12*c25+c15*c22-c15*c35+c20*c30-c30*c7;
t10 = -c11*c32-c12*c31+c16*c27+c17*c26;
t11 = -c11*c31+c16*c26;
t12 = -c11*c33-c13*c31-c16*c22+c16*c28-c17*c21+c18*c26+c31*c7+c32*c6;
t13 = -c11*c34-c14*c31+c16*c29+c19*c26;
t14 = c11*c22-c11*c35+c12*c21-c13*c34-c14*c33-c15*c31+c16*c30-c17*c24+c18*c29-c19*c22+c19*c28+c20*c26-c26*c7-c27*c6+c32*c9+c34*c7;
t15 = c14*c24-c29*c9;
t16 = c10*c33+c13*c23-c18*c25-c20*c23-c28*c8+c35*c8;
t17 = -c18*c23+c33*c8;
t18 = c12*c22-c12*c35-c15*c32+c17*c30+c20*c27-c27*c7;
t19 = -c12*c32+c17*c27;
t20 = -c16*c21+c31*c6;
t21 = -c16*c23-c18*c21+c31*c8+c33*c6;
t22 = c11*c21-c16*c24-c19*c21-c26*c6+c31*c9+c34*c6;
t23 = c10*c31+c11*c23+c13*c21-c16*c25-c18*c24-c19*c23-c20*c21-c26*c8-c28*c6+c33*c9+c34*c8+c35*c6;
t24 = c11*c24+c14*c21-c19*c24-c26*c9-c29*c6+c34*c9;
t25 = c12*c24+c14*c22-c14*c35-c15*c34+c19*c30+c20*c29-c27*c9-c29*c7;
t26 = -c10*c26+c10*c34+c11*c25+c13*c24+c14*c23+c15*c21-c19*c25-c20*c24-c28*c9-c29*c8-c30*c6+c35*c9;
t27 = -c12*c33-c13*c32-c17*c22+c17*c28+c18*c27+c32*c7;
t28 = -c12*c34-c14*c32+c17*c29+c19*c27;
t29 = -c13*c33-c17*c23-c18*c22+c18*c28+c32*c8+c33*c7;
t30 = c10*c32+c12*c23+c13*c22-c13*c35-c15*c33-c17*c25+c18*c30-c20*c22+c20*c28-c27*c8-c28*c7+c35*c7;
t31 = -c14*c34+c19*c29;
t32 = -c15*c20-c30*c35;
t33 = -c13*c20-c15^2-c15*c18+c20^2-c28*c35-c30^2-c30*c33+c35^2;
t34 = -c14*c20-c15*c19-c29*c35-c30*c34;
t35 = -c10*c15-c12*c20-c15*c17-c25*c30-c27*c35-c30*c32;
t36 = c16*c6+c21*c31;
t37 = c11*c17+c12*c16+c16*c8+c18*c6+c21*c33+c23*c31+c26*c32+c27*c31;
t38 = -c11*c7-c12*c6+c17*c9+c19*c7-c21*c27-c22*c26+c22*c34+c24*c32;
t39 = -c11*c6+c16*c9+c19*c6-c21*c26+c21*c34+c24*c31;
t40 = c10*c16-2*c11*c12-c11*c8+c12*c19-c13*c6+c14*c17+2*c16*c17+c18*c9+c19*c8+c20*c6-c21*c28+c21*c35-c23*c26+c23*c34+c24*c33+c25*c31-2*c26*c27+c27*c34+c29*c32+2*c31*c32;
t41 = -c11*c9-c14*c6+c19*c9-c21*c29-c24*c26+c24*c34;
t42 = -c10*c11+c10*c19-c11*c17-2*c12*c14-c12*c16-c13*c9-c14*c8-c15*c6+2*c17*c19+c20*c9-c21*c30-c23*c29-c24*c28+c24*c35-c25*c26+c25*c34-c26*c32-2*c27*c29-c27*c31+2*c32*c34;
t43 = c16*c7+c17*c6+c21*c32+c22*c31;
t44 = c13*c18+c28*c33;
t45 = -c14*c19-c29*c34;
t46 = -2*c13*c15-c13*c18+c15*c20+2*c18*c20-2*c28*c30-c28*c33+c30*c35+2*c33*c35;
t47 = -c13^2+c13*c20+c15*c18+c18^2-c28^2+c28*c35+c30*c33+c33^2;
t48 = -c10*c12-c12*c17-c15*c7-c22*c30-c25*c27-c27*c32;
t49 = -c12*c7-c22*c27;
t50 = c11*c16+c26*c31;
t51 = c11*c18+c13*c16+c26*c33+c28*c31;
t52 = c17*c7+c22*c32;
t53 = c12*c17+c17*c8+c18*c7+c22*c33+c23*c32+c27*c32;
t54 = c12*c18+c13*c17+c18*c8+c23*c33+c27*c33+c28*c32;
t55 = -c11^2+c11*c19+c14*c16+c16^2-c26^2+c26*c34+c29*c31+c31^2;
t56 = -2*c11*c13+c11*c20+c13*c19+c14*c18+c15*c16+2*c16*c18-2*c26*c28+c26*c35+c28*c34+c29*c33+c30*c31+2*c31*c33;
t57 = -2*c11*c14-c11*c16+c14*c19+2*c16*c19-2*c26*c29-c26*c31+c29*c34+2*c31*c34;
t58 = -2*c11*c15-c11*c18-2*c13*c14-c13*c16+c14*c20+c15*c19+2*c16*c20+2*c18*c19-2*c26*c30-c26*c33-2*c28*c29-c28*c31+c29*c35+c30*c34+2*c31*c35+2*c33*c34;
t59 = -c11*c19-c14^2-c14*c16+c19^2-c26*c34-c29^2-c29*c31+c34^2;
t60 = -c10*c14-c12*c19-c14*c17-c15*c9-c24*c30-c25*c29-c27*c34-c29*c32;
t61 = -c11*c20-c13*c19-2*c14*c15-c14*c18-c15*c16+2*c19*c20-c26*c35-c28*c34-2*c29*c30-c29*c33-c30*c31+2*c34*c35;
t62 = c10*c17-c12^2-c12*c8-c13*c7+c17^2+c20*c7-c22*c28+c22*c35-c23*c27+c25*c32-c27^2+c32^2;
t63 = -c12*c9-c14*c7-c22*c29-c24*c27;
t64 = c10*c18-2*c12*c13+c12*c20-c13*c8+c15*c17+2*c17*c18+c20*c8-c23*c28+c23*c35+c25*c33-2*c27*c28+c27*c35+c30*c32+2*c32*c33;
t65 = -c10*c13+c10*c20-2*c12*c15-c12*c18-c13*c17-c15*c8+2*c17*c20-c23*c30-c25*c28+c25*c35-2*c27*c30-c27*c33-c28*c32+2*c32*c35;
t66 = -c14*c9-c24*c29;
t67 = -c10*c20-c25*c35;
t68 = -c10*c15-c10*c18-c20*c8-c23*c35-c25*c30-c25*c33;
t69 = -c10*c19-c20*c9-c24*c35-c25*c34;
t70 = -c10^2-c10*c17+c20^2-c20*c7-c22*c35-c25^2-c25*c32+c35^2;
t71 = c16*c9+c19*c6+c21*c34+c24*c31;
t72 = c10*c16+c12*c19+c14*c17+2*c16*c17+c18*c9+c19*c8+c20*c6-2*c21*c22+c21*c35+c23*c34+c24*c33+c25*c31+c27*c34+c29*c32+2*c31*c32-2*c6*c7;
t73 = c11*c19+c14*c16+c16^2-c21^2+c26*c34+c29*c31+c31^2-c6^2;
t74 = c11*c20-c11*c7-c12*c6+c13*c19+c14*c18+c15*c16+2*c16*c18-2*c21*c23-c21*c27-c22*c26+c26*c35+c28*c34+c29*c33+c30*c31+2*c31*c33-2*c6*c8;
t75 = c14*c19+2*c16*c19-2*c21*c24+c29*c34+2*c31*c34-2*c6*c9;
t76 = -2*c10*c6-c12*c9+c14*c20-c14*c7+c15*c19+2*c16*c20-c16*c7-c17*c6+2*c18*c19-2*c21*c25-c21*c32-c22*c29-c22*c31-2*c23*c24-c24*c27+c29*c35+c30*c34+2*c31*c35+2*c33*c34-2*c8*c9;
t77 = -c19*c9-c24*c34;
t78 = -c10*c13-c15*c8-c18*c8-c23*c30-c23*c33-c25*c28;
t79 = -c13*c8-c23*c28;
t80 = c10*c20-2*c10*c7+2*c17*c20-c17*c7-2*c22*c25-c22*c32+c25*c35+2*c32*c35;
t81 = c10*c17+c17^2+c20*c7-c22^2+c22*c35+c25*c32+c32^2-c7^2;
t82 = c17*c9+c19*c7+c22*c34+c24*c32;
t83 = c19*c9+c24*c34;
t84 = -c11*c6-c21*c26;
t85 = -c11*c8-c13*c6-c21*c28-c23*c26;
t86 = -c11*c9-c14*c6-c16*c6-c21*c29-c21*c31-c24*c26;
t87 = -c10*c11-c13*c9-c14*c8-c15*c6-c16*c8-c18*c6-c21*c30-c21*c33-c23*c29-c23*c31-c24*c28-c25*c26;
t88 = -c14*c9-c16*c9-c19*c6-c21*c34-c24*c29-c24*c31;
t89 = -2*c10*c9-c17*c9+2*c19*c20-c19*c7-c22*c34-2*c24*c25-c24*c32+2*c34*c35;
t90 = -c10*c14-c10*c16-c15*c9-c18*c9-c19*c8-c20*c6-c21*c35-c23*c34-c24*c30-c24*c33-c25*c29-c25*c31;
t91 = c10*c18+c12*c20-c12*c7+c15*c17+2*c17*c18+c20*c8-2*c22*c23-c22*c27+c23*c35+c25*c33+c27*c35+c30*c32+2*c32*c33-2*c7*c8;
t92 = c10*c19+2*c17*c19+c20*c9-2*c22*c24+c24*c35+c25*c34+2*c32*c34-2*c7*c9;
t93 = -c12*c8+c13*c20-c13*c7+c15*c18+c18^2-c22*c28-c23^2-c23*c27+c28*c35+c30*c33+c33^2-c8^2;
t94 = -c10*c12-2*c10*c8+c15*c20-c15*c7-c17*c8+2*c18*c20-c18*c7-c22*c30-c22*c33-2*c23*c25-c23*c32-c25*c27+c30*c35+2*c33*c35;
t95 = c19^2-c24^2+c34^2-c9^2;
t96 = c15*c20+c30*c35;
t97 = c12*c7+c22*c27;
t98 = c11*c7+c12*c6+c21*c27+c22*c26;
t99 = c11*c6+c21*c26;
t100 = c11*c9+c14*c6+c21*c29+c24*c26;
t101 = c11^2-c21^2+c26^2-c6^2;
t102 = 2*c11*c13-c11*c7-c12*c6-2*c21*c23-c21*c27-c22*c26+2*c26*c28-2*c6*c8;
t103 = c10*c11+c11*c17+2*c12*c14+c12*c16+c13*c9+c14*c8+c15*c6+c21*c30-2*c22*c24+c23*c29+c24*c28+c25*c26+c26*c32+2*c27*c29+c27*c31-2*c7*c9;
t104 = 2*c11*c14+c11*c16-2*c21*c24+2*c26*c29+c26*c31-2*c6*c9;
t105 = -2*c10*c6+2*c11*c15+c11*c18-c12*c9+2*c13*c14+c13*c16-c14*c7-c16*c7-c17*c6-2*c21*c25-c21*c32-c22*c29-c22*c31-2*c23*c24-c24*c27+2*c26*c30+c26*c33+2*c28*c29+c28*c31-2*c8*c9;
t106 = c11*c19+c14^2+c14*c16-c24^2+c26*c34+c29^2+c29*c31-c9^2;
t107 = -2*c10*c9+c11*c20+c13*c19+2*c14*c15+c14*c18+c15*c16-c17*c9-c19*c7-c22*c34-2*c24*c25-c24*c32+c26*c35+c28*c34+2*c29*c30+c29*c33+c30*c31;
t108 = 2*c11*c12+c11*c8+c13*c6-2*c21*c22+c21*c28+c23*c26+2*c26*c27-2*c6*c7;
t109 = c10*c15+c12*c20+c15*c17+c25*c30+c27*c35+c30*c32;
t110 = c10*c12+c12*c17+c15*c7+c22*c30+c25*c27+c27*c32;
t111 = c12*c9+c14*c7+c22*c29+c24*c27;
t112 = c14*c9+c24*c29;
t113 = c12^2+c12*c8+c13*c7-c22^2+c22*c28+c23*c27+c27^2-c7^2;
t114 = 2*c12*c13-c12*c7+c13*c8-2*c22*c23-c22*c27+c23*c28+2*c27*c28-2*c7*c8;
t115 = -c12*c8+c13^2-c13*c7-c22*c28-c23^2-c23*c27+c28^2-c8^2;
t116 = c14*c20+c15*c19+c29*c35+c30*c34;
t117 = c10*c13-2*c10*c7+2*c12*c15+c12*c18+c13*c17+c15*c8-c17*c7-2*c22*c25-c22*c32+c23*c30+c25*c28+2*c27*c30+c27*c33+c28*c32;
t118 = c10*c14+c12*c19+c14*c17+c15*c9+c24*c30+c25*c29+c27*c34+c29*c32;
t119 = -c10*c12-2*c10*c8+2*c13*c15+c13*c18-c15*c7-c17*c8-c18*c7-c22*c30-c22*c33-2*c23*c25-c23*c32-c25*c27+2*c28*c30+c28*c33;
t120 = -c10^2-c10*c17+c13*c20+c15^2+c15*c18-c20*c7-c22*c35-c25^2-c25*c32+c28*c35+c30^2+c30*c33;
t121 = c14*c19+c29*c34;
M = zeros(6,42);
M(1,35) = -1;
M(1,36) = t2;
M(1,38) = t1;
M(1,39) = t3;
M(1,40) = t4;
M(1,42) = t5;
M(2,9) = t11;
M(2,10) = t20;
M(2,17) = t10;
M(2,18) = t12;
M(2,19) = t21;
M(2,21) = t13;
M(2,22) = t22;
M(2,23) = t19;
M(2,24) = t27;
M(2,25) = t29;
M(2,26) = t17;
M(2,27) = t28;
M(2,28) = t14;
M(2,29) = t23;
M(2,30) = t31;
M(2,31) = t24;
M(2,32) = t18;
M(2,33) = t30;
M(2,34) = t16;
M(2,35) = t25;
M(2,36) = t26;
M(2,37) = t15;
M(2,38) = t9;
M(2,39) = t7;
M(2,40) = t8;
M(2,42) = t6;
M(3,2) = t36;
M(3,3) = t50;
M(3,5) = t43;
M(3,6) = t37;
M(3,7) = t51;
M(3,9) = t39;
M(3,10) = t55;
M(3,12) = t52;
M(3,13) = t53;
M(3,14) = t54;
M(3,15) = t44;
M(3,17) = t38;
M(3,18) = t40;
M(3,19) = t56;
M(3,21) = t41;
M(3,22) = t57;
M(3,23) = t49;
M(3,24) = t62;
M(3,25) = t64;
M(3,26) = t47;
M(3,27) = t63;
M(3,28) = t42;
M(3,29) = t58;
M(3,30) = t66;
M(3,31) = t59;
M(3,32) = t48;
M(3,33) = t65;
M(3,34) = t46;
M(3,35) = t60;
M(3,36) = t61;
M(3,37) = t45;
M(3,38) = t35;
M(3,39) = t33;
M(3,40) = t34;
M(3,42) = t32;
M(4,1) = t36;
M(4,2) = t50;
M(4,4) = t43;
M(4,5) = t37;
M(4,6) = t51;
M(4,8) = t71;
M(4,9) = t73;
M(4,10) = t84;
M(4,11) = t52;
M(4,12) = t53;
M(4,13) = t54;
M(4,14) = t44;
M(4,16) = t82;
M(4,17) = t72;
M(4,18) = t74;
M(4,19) = t85;
M(4,20) = t83;
M(4,21) = t75;
M(4,22) = t86;
M(4,23) = t81;
M(4,24) = t91;
M(4,25) = t93;
M(4,26) = t79;
M(4,27) = t92;
M(4,28) = t76;
M(4,29) = t87;
M(4,30) = t95;
M(4,31) = t88;
M(4,32) = t80;
M(4,33) = t94;
M(4,34) = t78;
M(4,35) = t89;
M(4,36) = t90;
M(4,37) = t77;
M(4,38) = t70;
M(4,39) = t68;
M(4,40) = t69;
M(4,42) = t67;
M(5,1) = t99;
M(5,2) = t101;
M(5,3) = t84;
M(5,4) = t98;
M(5,5) = t108;
M(5,6) = t102;
M(5,7) = t85;
M(5,8) = t100;
M(5,9) = t104;
M(5,10) = t86;
M(5,11) = t97;
M(5,12) = t113;
M(5,13) = t114;
M(5,14) = t115;
M(5,15) = t79;
M(5,16) = t111;
M(5,17) = t103;
M(5,18) = t105;
M(5,19) = t87;
M(5,20) = t112;
M(5,21) = t106;
M(5,22) = t88;
M(5,23) = t110;
M(5,24) = t117;
M(5,25) = t119;
M(5,26) = t78;
M(5,27) = t118;
M(5,28) = t107;
M(5,29) = t90;
M(5,30) = t121;
M(5,31) = t77;
M(5,32) = t109;
M(5,33) = t120;
M(5,34) = t68;
M(5,35) = t116;
M(5,36) = t69;
M(5,38) = t96;
M(5,39) = t67;
M(6,40) = 1;
M(6,41) = -1;
M = [rref(M(1:0,:)); M(1:end,:)];
ncinds = [4, 5, 9, 10, 11, 15, 17, 22, 23, 27, 28, 29, 33, 34, 35, 39, 41, 46, 47, 50, 51, 52, 53, 56, 57, 58, 59, 64, 65, 69, 70, 71, 75, 76, 77, 81, 82, 83, 87, 89, 94, 95, 98, 99, 100, 101, 104, 105, 106, 107, 110, 111, 112, 113, 118, 119, 122, 123, 124, 125, 128, 129, 130, 131, 134, 135, 136, 137, 140, 141, 142, 143, 146, 147, 148, 149, 152, 153, 154, 155, 158, 159, 160, 161, 164, 165, 166, 167, 170, 171, 172, 173, 176, 177, 178, 179, 182, 183, 184, 185, 188, 189, 190, 191, 194, 195, 196, 197, 200, 201, 202, 203, 206, 207, 208, 209, 211, 212, 213, 214, 215, 218, 219, 220, 223, 224, 225, 226, 227, 229, 230, 231, 232, 233, 235, 236, 237, 238, 247, 248, 249, 250];
for ncind = 1:142
    eval(strjoin({'nc', num2str(ncind), ' = ', 'M(ncinds(ncind));'}, '') )
end
ArrayOfCsNames = '';
sizeOfC = 41;
noOfVars = 4;
sparseBasis = [1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 1 1 1 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3; 0 0 0 0 1 1 1 2 2 2 3 3 0 0 0 1 1 2 2 0 1 4 0 0 1 1 2 3 3 0 0 0 0 0 1 1 1 1 2 2 2; 0 2 3 4 1 2 3 0 1 2 0 1 0 1 2 0 2 1 2 1 0 0 3 4 1 3 0 0 1 0 1 2 3 4 0 1 2 3 0 1 2;];
solFromEigenVectors = [-1 -0.1e1 / 0.2e1 -0.1e1 / 0.2e1; 0 -0.1e1 / 0.2e1 0.1e1 / 0.2e1; 0 0 0; 0 0 0; 0 1 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 1 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0;];
Cs = zeros(41,82);
Cs(1,1) = nc139;
Cs(1,13) = nc135;
Cs(1,14) = nc117;
Cs(1,16) = -1;
Cs(1,20) = nc130;
Cs(1,21) = nc125;
Cs(2,2) = nc130;
Cs(2,5) = nc125;
Cs(2,14) = nc135;
Cs(2,15) = nc117;
Cs(2,20) = nc139;
Cs(2,25) = -1;
Cs(3,2) = nc139;
Cs(3,3) = nc130;
Cs(3,6) = nc125;
Cs(3,15) = nc135;
Cs(3,17) = -1;
Cs(3,23) = nc117;
Cs(4,3) = nc139;
Cs(4,4) = nc130;
Cs(4,7) = nc125;
Cs(4,23) = nc135;
Cs(4,24) = nc117;
Cs(4,26) = -1;
Cs(5,5) = nc130;
Cs(5,8) = nc125;
Cs(5,16) = nc135;
Cs(5,21) = nc139;
Cs(5,25) = nc117;
Cs(5,27) = -1;
Cs(6,5) = nc139;
Cs(6,6) = nc130;
Cs(6,9) = nc125;
Cs(6,17) = nc117;
Cs(6,18) = -1;
Cs(6,25) = nc135;
Cs(7,6) = nc139;
Cs(7,7) = nc130;
Cs(7,10) = nc125;
Cs(7,17) = nc135;
Cs(7,19) = -1;
Cs(7,26) = nc117;
Cs(8,8) = nc139;
Cs(8,9) = nc130;
Cs(8,11) = nc125;
Cs(8,18) = nc117;
Cs(8,27) = nc135;
Cs(8,28) = -1;
Cs(9,9) = nc139;
Cs(9,10) = nc130;
Cs(9,12) = nc125;
Cs(9,18) = nc135;
Cs(9,19) = nc117;
Cs(9,29) = -1;
Cs(10,13) = nc139;
Cs(10,14) = nc130;
Cs(10,16) = nc125;
Cs(10,30) = nc135;
Cs(10,31) = nc117;
Cs(10,35) = -1;
Cs(11,14) = nc139;
Cs(11,15) = nc130;
Cs(11,25) = nc125;
Cs(11,31) = nc135;
Cs(11,32) = nc117;
Cs(11,36) = -1;
Cs(12,15) = nc139;
Cs(12,17) = nc125;
Cs(12,23) = nc130;
Cs(12,32) = nc135;
Cs(12,33) = nc117;
Cs(12,37) = -1;
Cs(13,23) = nc139;
Cs(13,24) = nc130;
Cs(13,26) = nc125;
Cs(13,33) = nc135;
Cs(13,34) = nc117;
Cs(13,38) = -1;
Cs(14,16) = nc139;
Cs(14,25) = nc130;
Cs(14,27) = nc125;
Cs(14,35) = nc135;
Cs(14,36) = nc117;
Cs(14,39) = -1;
Cs(15,17) = nc130;
Cs(15,18) = nc125;
Cs(15,25) = nc139;
Cs(15,36) = nc135;
Cs(15,37) = nc117;
Cs(15,40) = -1;
Cs(16,17) = nc139;
Cs(16,19) = nc125;
Cs(16,26) = nc130;
Cs(16,37) = nc135;
Cs(16,38) = nc117;
Cs(16,41) = -1;
Cs(17,1) = nc140;
Cs(17,2) = nc109;
Cs(17,3) = nc77;
Cs(17,5) = nc105;
Cs(17,6) = nc73;
Cs(17,8) = nc101;
Cs(17,9) = nc69;
Cs(17,11) = nc65;
Cs(17,13) = nc136;
Cs(17,14) = nc118;
Cs(17,15) = nc89;
Cs(17,16) = nc113;
Cs(17,17) = nc47;
Cs(17,18) = nc43;
Cs(17,20) = nc131;
Cs(17,21) = nc126;
Cs(17,23) = nc51;
Cs(17,25) = nc85;
Cs(17,27) = nc81;
Cs(17,30) = nc122;
Cs(17,31) = nc97;
Cs(17,32) = nc61;
Cs(17,33) = nc24;
Cs(17,35) = nc93;
Cs(17,36) = nc57;
Cs(17,37) = nc20;
Cs(18,2) = nc131;
Cs(18,3) = nc109;
Cs(18,4) = nc77;
Cs(18,5) = nc126;
Cs(18,6) = nc105;
Cs(18,7) = nc73;
Cs(18,9) = nc101;
Cs(18,10) = nc69;
Cs(18,12) = nc65;
Cs(18,14) = nc136;
Cs(18,15) = nc118;
Cs(18,17) = nc85;
Cs(18,18) = nc81;
Cs(18,19) = nc43;
Cs(18,20) = nc140;
Cs(18,23) = nc89;
Cs(18,24) = nc51;
Cs(18,25) = nc113;
Cs(18,26) = nc47;
Cs(18,31) = nc122;
Cs(18,32) = nc97;
Cs(18,33) = nc61;
Cs(18,34) = nc24;
Cs(18,36) = nc93;
Cs(18,37) = nc57;
Cs(18,38) = nc20;
Cs(19,5) = nc131;
Cs(19,6) = nc109;
Cs(19,7) = nc77;
Cs(19,8) = nc126;
Cs(19,9) = nc105;
Cs(19,10) = nc73;
Cs(19,11) = nc101;
Cs(19,12) = nc69;
Cs(19,16) = nc136;
Cs(19,17) = nc89;
Cs(19,18) = nc85;
Cs(19,19) = nc47;
Cs(19,21) = nc140;
Cs(19,22) = nc65;
Cs(19,25) = nc118;
Cs(19,26) = nc51;
Cs(19,27) = nc113;
Cs(19,28) = nc81;
Cs(19,29) = nc43;
Cs(19,35) = nc122;
Cs(19,36) = nc97;
Cs(19,37) = nc61;
Cs(19,38) = nc24;
Cs(19,39) = nc93;
Cs(19,40) = nc57;
Cs(19,41) = nc20;
Cs(20,1) = nc141;
Cs(20,2) = nc110;
Cs(20,3) = nc78;
Cs(20,4) = nc39;
Cs(20,5) = nc106;
Cs(20,6) = nc74;
Cs(20,7) = nc36;
Cs(20,8) = nc102;
Cs(20,9) = nc70;
Cs(20,10) = nc33;
Cs(20,11) = nc66;
Cs(20,12) = nc30;
Cs(20,13) = nc137;
Cs(20,14) = nc119;
Cs(20,15) = nc90;
Cs(20,16) = nc114;
Cs(20,17) = nc48;
Cs(20,18) = nc44;
Cs(20,19) = nc10;
Cs(20,20) = nc132;
Cs(20,21) = nc127;
Cs(20,23) = nc52;
Cs(20,24) = nc16;
Cs(20,25) = nc86;
Cs(20,26) = nc13;
Cs(20,27) = nc82;
Cs(20,30) = nc123;
Cs(20,31) = nc98;
Cs(20,32) = nc62;
Cs(20,33) = nc25;
Cs(20,34) = nc6;
Cs(20,35) = nc94;
Cs(20,36) = nc58;
Cs(20,37) = nc21;
Cs(20,38) = nc3;
Cs(21,1) = nc142;
Cs(21,2) = nc111;
Cs(21,3) = nc79;
Cs(21,5) = nc107;
Cs(21,6) = nc75;
Cs(21,7) = nc37;
Cs(21,8) = nc103;
Cs(21,9) = nc71;
Cs(21,10) = nc34;
Cs(21,11) = nc67;
Cs(21,12) = nc31;
Cs(21,13) = nc138;
Cs(21,14) = nc120;
Cs(21,15) = nc91;
Cs(21,16) = nc115;
Cs(21,17) = nc49;
Cs(21,18) = nc45;
Cs(21,19) = nc11;
Cs(21,20) = nc133;
Cs(21,21) = nc128;
Cs(21,22) = nc28;
Cs(21,23) = nc53;
Cs(21,25) = nc87;
Cs(21,26) = nc14;
Cs(21,27) = nc83;
Cs(21,28) = nc41;
Cs(21,29) = nc8;
Cs(21,30) = nc124;
Cs(21,31) = nc99;
Cs(21,32) = nc63;
Cs(21,33) = nc26;
Cs(21,35) = nc95;
Cs(21,36) = nc59;
Cs(21,37) = nc22;
Cs(21,38) = nc4;
Cs(21,39) = nc55;
Cs(21,40) = nc18;
Cs(21,41) = nc1;
Cs(22,2) = nc112;
Cs(22,3) = nc80;
Cs(22,4) = nc40;
Cs(22,5) = nc108;
Cs(22,6) = nc76;
Cs(22,7) = nc38;
Cs(22,8) = nc104;
Cs(22,9) = nc72;
Cs(22,10) = nc35;
Cs(22,11) = nc68;
Cs(22,12) = nc32;
Cs(22,14) = nc121;
Cs(22,15) = nc92;
Cs(22,16) = nc116;
Cs(22,17) = nc50;
Cs(22,18) = nc46;
Cs(22,19) = nc12;
Cs(22,20) = nc134;
Cs(22,21) = nc129;
Cs(22,22) = nc29;
Cs(22,23) = nc54;
Cs(22,24) = nc17;
Cs(22,25) = nc88;
Cs(22,26) = nc15;
Cs(22,27) = nc84;
Cs(22,28) = nc42;
Cs(22,29) = nc9;
Cs(22,31) = nc100;
Cs(22,32) = nc64;
Cs(22,33) = nc27;
Cs(22,34) = nc7;
Cs(22,35) = nc96;
Cs(22,36) = nc60;
Cs(22,37) = nc23;
Cs(22,38) = nc5;
Cs(22,39) = nc56;
Cs(22,40) = nc19;
Cs(22,41) = nc2;
Cs(23,13) = 1;
Cs(23,42) = -1;
Cs(24,15) = 1;
Cs(24,43) = -1;
Cs(25,23) = 1;
Cs(25,44) = -1;
Cs(26,24) = 1;
Cs(26,45) = -1;
Cs(27,25) = 1;
Cs(27,46) = -1;
Cs(28,17) = 1;
Cs(28,47) = -1;
Cs(29,26) = 1;
Cs(29,48) = -1;
Cs(30,27) = 1;
Cs(30,49) = -1;
Cs(31,18) = 1;
Cs(31,50) = -1;
Cs(32,19) = 1;
Cs(32,51) = -1;
Cs(33,28) = 1;
Cs(33,52) = -1;
Cs(34,29) = 1;
Cs(34,53) = -1;
Cs(35,30) = 1;
Cs(35,54) = -1;
Cs(36,31) = 1;
Cs(36,55) = -1;
Cs(37,32) = 1;
Cs(37,56) = -1;
Cs(38,35) = 1;
Cs(38,57) = -1;
Cs(39,37) = 1;
Cs(39,58) = -1;
Cs(40,40) = 1;
Cs(40,59) = -1;
Cs(41,41) = 1;
Cs(41,60) = -1;
indicesToRemove = [];
indicesToSkip = fliplr(find(sum(transpose(solFromEigenVectors))~=0));
sizeoffinalres =19
rowsToRemove = [];
colsToRemove = [];
rowsToRemove = rowsToRemove(1:max(find(rowsToRemove)), :);
colsToRemove = colsToRemove(1:max(find(colsToRemove)), :);
% if(size(rowsToRemove,1) > 0)
if 1 == 2
    indicesToRemove = [rowsToRemove,colsToRemove];
    
    noOfCoeffMatrices = size(Cs,2)/size(Cs,1);
    
    sizeOfReducedCs = sizeOfC - size(rowsToRemove,1);
    sizeOfM22 = size(rowsToRemove,1);
    numberOfCoeffMs = size(Cs,2)/size(Cs,1);
    maxReducedCoeffMs = 2 * numberOfCoeffMs - 1;
    colsToRemove = colsToRemove + ([1:numberOfCoeffMs] - 1) * sizeOfC;
    colsToRemove = colsToRemove(:);
    
    % Estimating M's
    M11 = Cs;
    M11(rowsToRemove,:) =[];
    M11(:,colsToRemove) =[];
    M22 = Cs(rowsToRemove, colsToRemove);
    M12 = Cs(:,colsToRemove);
    M12(rowsToRemove,:) = [];
    M21 = Cs(rowsToRemove,:);
    M21(:,colsToRemove) = [];
    reducedCs = [M11,zeros(sizeOfReducedCs,sizeOfReducedCs*(maxReducedCoeffMs-numberOfCoeffMs))];
    M22 = M22(1:sizeOfM22, 1:sizeOfM22);
    invM22 = inv(M22);
    indxOfZeroCoeffMat = [];
    for i=1:numberOfCoeffMs
        for j=1:numberOfCoeffMs
            coeffExp = i+j-1;
            temp = reducedCs(:,sizeOfReducedCs*(coeffExp-1) + 1:sizeOfReducedCs*coeffExp);
            tempM12 = M12(:,sizeOfM22*(i-1) + 1:sizeOfM22*i);
            tempM21 = M21(:,sizeOfReducedCs*(j-1) + 1:sizeOfReducedCs*j);
            
            buff = temp - (tempM12*invM22) * tempM21;
            reducedCs(:,sizeOfReducedCs*(coeffExp-1) + 1:sizeOfReducedCs*coeffExp) = buff;
            if( norm(buff,'fro') ~= 0)
                indxOfZeroCoeffMat(coeffExp) = 1;
            else
                indxOfZeroCoeffMat(coeffExp) = 0;
            end
        end
    end
    
    allReducedCs = [];
    for i = 1:maxReducedCoeffMs
        buff = reducedCs(:,sizeOfReducedCs*(i-1) + 1:sizeOfReducedCs*i);
        if(indxOfZeroCoeffMat(i) == 1)
            allReducedCs = [allReducedCs, buff];
            
        end
    end
    numOfRedCoeffMat = max(find(indxOfZeroCoeffMat));
    solForm = solFromEigenVectors;
    solForm(indicesToRemove(:,2),:) = [];
    sparseBasis(:,indicesToRemove(:,2)) = [];
    allCss = mat2cell(allReducedCs, sizeOfReducedCs, ones(1,numOfRedCoeffMat) * sizeOfReducedCs);
    
else
    sizeOfReducedCs = sizeOfC;
    numOfRedCoeffMat = size(Cs,2)/size(Cs,1);
    allCss = mat2cell(Cs, sizeOfReducedCs, ones(1,numOfRedCoeffMat ) * sizeOfReducedCs);
    indicesToRemove = [];
end


if rank(allCss{end}) < rank(allCss{1})
    allCss = fliplr(allCss);
end

A = []; B = [];
for i = 1:numOfRedCoeffMat -2
    tempA = []; tempB = [];
    for j = 1:numOfRedCoeffMat - 1
        if j == i + 1
            tempA = [tempA eye(sizeOfReducedCs)];
        else
            tempA = [tempA zeros(sizeOfReducedCs)];
        end
        if j == i
            tempB = [tempB eye(sizeOfReducedCs)];
        else
            tempB = [tempB zeros(sizeOfReducedCs)];
        end
    end
    A = [A;tempA];
    B = [B;tempB];
end

tempA = []; tempB = [];
for j = 1:numOfRedCoeffMat - 1
    tempA = [tempA -allCss{j}];
    if j == numOfRedCoeffMat -1
        tempB = [tempB allCss{j+1}];
    else
        tempB = [tempB zeros(sizeOfReducedCs)];
    end
end
A = [A; tempA];
B = [B; tempB];



extendedbasis = [sparseBasis; zeros(1,sizeOfReducedCs)];
for i = 1:numOfRedCoeffMat-2
    extendedbasis = [extendedbasis, [sparseBasis;ones(1,sizeOfReducedCs)*i]];
end


%% Removal of 0 eigen values via row algebra operations.

perms = [];
zeigvalindx=[];
colpermutation = [];
rowpermutation = [];

% remove_0_eigvalues;
% remove_more_0_eigvalues;


%% Generation of the template to extract solutions to individual
% unhidden variables from the eigen vector which will be a vector
% of monomials, which is an extension of the sparse bases.

% sizeoffinalres = size(extendedbasis,2);
% sizeoffinalres = min(sizeoffinalres, size(extendedbasis,2));
tempextendedbasis = extendedbasis(:,1:sizeoffinalres);
% tempextendedbasis = extendedbasis(:,end-sizeoffinalres+1:end);    
disp(tempextendedbasis);
numerr = 1e-10;

C0 = B; C1 = A;

A1 = -C0(1:end-sizeoffinalres, 1:end-sizeoffinalres);
A2 = -C0(1:end-sizeoffinalres, end-sizeoffinalres+1:end);
B1 = C1(end-sizeoffinalres+1:end, 1:end-sizeoffinalres);
B2 = C1(end-sizeoffinalres+1:end, end-sizeoffinalres+1:end);
X = B2 - (B1/A1) * A2;

[X, nullspacesize, depdXcols, nzrows, depdCcols, indepdCcols, mdepdXcols, mindepdXcols, indepdCrows, rowstorem] = remove_more_0_eigvalues_by_deflation(X);
nullspacesize1=nullspacesize;
depdXcols1=depdXcols;
nzrows1=nzrows;
depdCcols1=depdCcols;
indepdCcols1=indepdCcols;
mdepdXcols1=mdepdXcols;
mindepdXcols1=mindepdXcols;
indepdCrows1=indepdCrows;
rowstorem1=rowstorem;

[X, nullspacesize, depdXcols, nzrows, depdCcols, indepdCcols, mdepdXcols, mindepdXcols, indepdCrows, rowstorem] = remove_more_0_eigvalues_by_deflation(X);
nullspacesize2=nullspacesize;
depdXcols2=depdXcols;
nzrows2=nzrows;
depdCcols2=depdCcols;
indepdCcols2=indepdCcols;
mdepdXcols2=mdepdXcols;
mindepdXcols2=mindepdXcols;
indepdCrows2=indepdCrows;
rowstorem2=rowstorem;

[X, nullspacesize, depdXcols, nzrows, depdCcols, indepdCcols, mdepdXcols, mindepdXcols, indepdCrows, rowstorem] = remove_more_0_eigvalues_by_deflation(X);
nullspacesize3=nullspacesize;
depdXcols3=depdXcols;
nzrows3=nzrows;
depdCcols3=depdCcols;
indepdCcols3=indepdCcols;
mdepdXcols3=mdepdXcols;
mindepdXcols3=mindepdXcols;
indepdCrows3=indepdCrows;
rowstorem3=rowstorem;
    

pairs = transpose(combnk(setdiff(1:size(tempextendedbasis,2), rowstorem),2));

ops = eye(size(tempextendedbasis,1));
solverTemplate = zeros(size(tempextendedbasis));
for i = 1:size(tempextendedbasis,1) -1
    op = ops(:,i);
    for pair = pairs
        tmp = tempextendedbasis(:,pair)* [1;-1];
        if norm(abs(tmp) - op) == 0
            solverTemplate(:,pair(1)) = solverTemplate(:,pair(1)) + sum(tmp) * op;
            solverTemplate(:,pair(2)) = solverTemplate(:,pair(2)) -sum(tmp) * op;
            break;
        end
    end
end
%     Removoing the extended part of the sparseBasis
extendedbasis = tempextendedbasis;
%%


solverTemplate(size(extendedbasis,1),:) = [];

disp(" Parasitic eigenvalues have been removed. And size of final solver is ");
disp(size(extendedbasis));
end
