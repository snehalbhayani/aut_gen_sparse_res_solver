function[PEPsolutions] = solve(data) 
hiddenvarnumber = 3;
c1 = data(1);c2 = data(2);c3 = data(3);c4 = data(4);c5 = data(5);c6 = data(6);c7 = data(7);c8 = data(8);c9 = data(9);c10 = data(10);c11 = data(11);c12 = data(12);c13 = data(13);c14 = data(14);c15 = data(15);c16 = data(16);c17 = data(17);c18 = data(18);c19 = data(19);c20 = data(20);c21 = data(21);c22 = data(22);c23 = data(23);c24 = data(24);c25 = data(25);c26 = data(26);c27 = data(27);c28 = data(28);c29 = data(29);c30 = data(30);c31 = data(31);c32 = data(32);c33 = data(33);c34 = data(34);c35 = data(35);
t1 = 4*c21;
t2 = 3*c22;
t3 = 3*c23;
t4 = 3*c11;
t5 = 2*c24;
t6 = 2*c25;
t7 = 2*c12;
t8 = 2*c26;
t9 = 2*c13;
t10 = 2*c5;
t11 = c27;
t12 = c28;
t13 = c14;
t14 = c29;
t15 = c15;
t16 = c6;
t17 = c30;
t18 = c16;
t19 = c7;
t20 = c2;
t21 = c22;
t22 = c25;
t23 = c12;
t24 = 3*c27;
t25 = 2*c28;
t26 = 2*c14;
t27 = 4*c31;
t28 = 3*c32;
t29 = 3*c17;
t30 = 2*c33;
t31 = 2*c18;
t32 = 2*c8;
t33 = c34;
t34 = c19;
t35 = c9;
t36 = c3;
t37 = c23;
t38 = c13;
t39 = 2*c29;
t40 = 3*c30;
t41 = 2*c16;
t42 = c32;
t43 = c18;
t44 = 3*c34;
t45 = 2*c19;
t46 = 4*c35;
t47 = 3*c20;
t48 = 2*c10;
t49 = c4;
M = zeros(4,21);
M(1,1) = t17;
M(1,2) = t8;
M(1,3) = t3;
M(1,4) = t1;
M(1,5) = t14;
M(1,6) = t6;
M(1,7) = t2;
M(1,8) = t12;
M(1,9) = t5;
M(1,10) = t11;
M(1,11) = t18;
M(1,12) = t9;
M(1,13) = t4;
M(1,14) = t15;
M(1,15) = t7;
M(1,16) = t13;
M(1,17) = t19;
M(1,18) = t10;
M(1,19) = t16;
M(1,21) = t20;
M(2,1) = t33;
M(2,2) = t14;
M(2,3) = t22;
M(2,4) = t21;
M(2,5) = t30;
M(2,6) = t25;
M(2,7) = t5;
M(2,8) = t28;
M(2,9) = t24;
M(2,10) = t27;
M(2,11) = t34;
M(2,12) = t15;
M(2,13) = t23;
M(2,14) = t31;
M(2,15) = t26;
M(2,16) = t29;
M(2,17) = t35;
M(2,18) = t16;
M(2,19) = t32;
M(2,21) = t36;
M(3,1) = t46;
M(3,2) = t40;
M(3,3) = t8;
M(3,4) = t37;
M(3,5) = t44;
M(3,6) = t39;
M(3,7) = t22;
M(3,8) = t30;
M(3,9) = t12;
M(3,10) = t42;
M(3,11) = t47;
M(3,12) = t41;
M(3,13) = t38;
M(3,14) = t45;
M(3,15) = t15;
M(3,16) = t43;
M(3,17) = t48;
M(3,18) = t19;
M(3,19) = t35;
M(3,21) = t49;
M(4,17) = 1;
M(4,20) = -1;
M = [rref(M(1:0,:)); M(1:end,:)];
nc1 = M(1);nc2 = M(2);nc3 = M(3);nc4 = M(5);nc5 = M(6);nc6 = M(7);nc7 = M(9);nc8 = M(10);nc9 = M(11);nc10 = M(13);nc11 = M(14);nc12 = M(15);nc13 = M(17);nc14 = M(18);nc15 = M(19);nc16 = M(21);nc17 = M(22);nc18 = M(23);nc19 = M(25);nc20 = M(26);nc21 = M(27);nc22 = M(29);nc23 = M(30);nc24 = M(31);nc25 = M(33);nc26 = M(34);nc27 = M(35);nc28 = M(37);nc29 = M(38);nc30 = M(39);nc31 = M(41);nc32 = M(42);nc33 = M(43);nc34 = M(45);nc35 = M(46);nc36 = M(47);nc37 = M(49);nc38 = M(50);nc39 = M(51);nc40 = M(53);nc41 = M(54);nc42 = M(55);nc43 = M(57);nc44 = M(58);nc45 = M(59);nc46 = M(61);nc47 = M(62);nc48 = M(63);nc49 = M(65);nc50 = M(66);nc51 = M(67);nc52 = M(69);nc53 = M(70);nc54 = M(71);nc55 = M(73);nc56 = M(74);nc57 = M(75);nc58 = M(81);nc59 = M(82);nc60 = M(83);
Cs = zeros(114,228);
Cs(1,1) = nc58;
Cs(1,2) = nc55;
Cs(1,3) = nc46;
Cs(1,4) = nc28;
Cs(1,7) = nc52;
Cs(1,8) = nc43;
Cs(1,9) = nc25;
Cs(1,12) = nc37;
Cs(1,13) = nc19;
Cs(1,17) = nc10;
Cs(1,27) = nc49;
Cs(1,34) = nc40;
Cs(1,35) = nc22;
Cs(1,40) = nc34;
Cs(1,41) = nc16;
Cs(1,46) = nc7;
Cs(1,61) = nc31;
Cs(1,62) = nc13;
Cs(1,67) = nc4;
Cs(1,82) = nc1;
Cs(2,2) = nc58;
Cs(2,3) = nc55;
Cs(2,4) = nc46;
Cs(2,5) = nc28;
Cs(2,8) = nc52;
Cs(2,9) = nc43;
Cs(2,10) = nc25;
Cs(2,13) = nc37;
Cs(2,14) = nc19;
Cs(2,18) = nc10;
Cs(2,34) = nc49;
Cs(2,35) = nc40;
Cs(2,36) = nc22;
Cs(2,41) = nc34;
Cs(2,42) = nc16;
Cs(2,47) = nc7;
Cs(2,62) = nc31;
Cs(2,63) = nc13;
Cs(2,68) = nc4;
Cs(2,83) = nc1;
Cs(3,3) = nc58;
Cs(3,4) = nc55;
Cs(3,5) = nc46;
Cs(3,6) = nc28;
Cs(3,9) = nc52;
Cs(3,10) = nc43;
Cs(3,11) = nc25;
Cs(3,14) = nc37;
Cs(3,15) = nc19;
Cs(3,19) = nc10;
Cs(3,35) = nc49;
Cs(3,36) = nc40;
Cs(3,37) = nc22;
Cs(3,42) = nc34;
Cs(3,43) = nc16;
Cs(3,48) = nc7;
Cs(3,63) = nc31;
Cs(3,64) = nc13;
Cs(3,69) = nc4;
Cs(3,84) = nc1;
Cs(4,7) = nc58;
Cs(4,8) = nc55;
Cs(4,9) = nc46;
Cs(4,10) = nc28;
Cs(4,12) = nc52;
Cs(4,13) = nc43;
Cs(4,14) = nc25;
Cs(4,17) = nc37;
Cs(4,18) = nc19;
Cs(4,21) = nc10;
Cs(4,40) = nc49;
Cs(4,41) = nc40;
Cs(4,42) = nc22;
Cs(4,46) = nc34;
Cs(4,47) = nc16;
Cs(4,51) = nc7;
Cs(4,67) = nc31;
Cs(4,68) = nc13;
Cs(4,72) = nc4;
Cs(4,87) = nc1;
Cs(5,8) = nc58;
Cs(5,9) = nc55;
Cs(5,10) = nc46;
Cs(5,11) = nc28;
Cs(5,13) = nc52;
Cs(5,14) = nc43;
Cs(5,15) = nc25;
Cs(5,18) = nc37;
Cs(5,19) = nc19;
Cs(5,22) = nc10;
Cs(5,41) = nc49;
Cs(5,42) = nc40;
Cs(5,43) = nc22;
Cs(5,47) = nc34;
Cs(5,48) = nc16;
Cs(5,52) = nc7;
Cs(5,68) = nc31;
Cs(5,69) = nc13;
Cs(5,73) = nc4;
Cs(5,88) = nc1;
Cs(6,12) = nc58;
Cs(6,13) = nc55;
Cs(6,14) = nc46;
Cs(6,15) = nc28;
Cs(6,17) = nc52;
Cs(6,18) = nc43;
Cs(6,19) = nc25;
Cs(6,21) = nc37;
Cs(6,22) = nc19;
Cs(6,24) = nc10;
Cs(6,46) = nc49;
Cs(6,47) = nc40;
Cs(6,48) = nc22;
Cs(6,51) = nc34;
Cs(6,52) = nc16;
Cs(6,55) = nc7;
Cs(6,72) = nc31;
Cs(6,73) = nc13;
Cs(6,76) = nc4;
Cs(6,91) = nc1;
Cs(7,13) = nc58;
Cs(7,14) = nc55;
Cs(7,15) = nc46;
Cs(7,16) = nc28;
Cs(7,18) = nc52;
Cs(7,19) = nc43;
Cs(7,20) = nc25;
Cs(7,22) = nc37;
Cs(7,23) = nc19;
Cs(7,25) = nc10;
Cs(7,47) = nc49;
Cs(7,48) = nc40;
Cs(7,49) = nc22;
Cs(7,52) = nc34;
Cs(7,53) = nc16;
Cs(7,56) = nc7;
Cs(7,73) = nc31;
Cs(7,74) = nc13;
Cs(7,77) = nc4;
Cs(7,92) = nc1;
Cs(8,14) = nc58;
Cs(8,15) = nc55;
Cs(8,16) = nc46;
Cs(8,19) = nc52;
Cs(8,20) = nc43;
Cs(8,23) = nc37;
Cs(8,28) = nc28;
Cs(8,29) = nc25;
Cs(8,30) = nc19;
Cs(8,31) = nc10;
Cs(8,48) = nc49;
Cs(8,49) = nc40;
Cs(8,50) = nc22;
Cs(8,53) = nc34;
Cs(8,54) = nc16;
Cs(8,57) = nc7;
Cs(8,74) = nc31;
Cs(8,75) = nc13;
Cs(8,78) = nc4;
Cs(8,93) = nc1;
Cs(9,17) = nc58;
Cs(9,18) = nc55;
Cs(9,19) = nc46;
Cs(9,20) = nc28;
Cs(9,21) = nc52;
Cs(9,22) = nc43;
Cs(9,23) = nc25;
Cs(9,24) = nc37;
Cs(9,25) = nc19;
Cs(9,26) = nc10;
Cs(9,51) = nc49;
Cs(9,52) = nc40;
Cs(9,53) = nc22;
Cs(9,55) = nc34;
Cs(9,56) = nc16;
Cs(9,58) = nc7;
Cs(9,76) = nc31;
Cs(9,77) = nc13;
Cs(9,79) = nc4;
Cs(9,94) = nc1;
Cs(10,18) = nc58;
Cs(10,19) = nc55;
Cs(10,20) = nc46;
Cs(10,22) = nc52;
Cs(10,23) = nc43;
Cs(10,25) = nc37;
Cs(10,29) = nc28;
Cs(10,30) = nc25;
Cs(10,31) = nc19;
Cs(10,32) = nc10;
Cs(10,52) = nc49;
Cs(10,53) = nc40;
Cs(10,54) = nc22;
Cs(10,56) = nc34;
Cs(10,57) = nc16;
Cs(10,59) = nc7;
Cs(10,77) = nc31;
Cs(10,78) = nc13;
Cs(10,80) = nc4;
Cs(10,95) = nc1;
Cs(11,21) = nc58;
Cs(11,22) = nc55;
Cs(11,23) = nc46;
Cs(11,24) = nc52;
Cs(11,25) = nc43;
Cs(11,26) = nc37;
Cs(11,30) = nc28;
Cs(11,31) = nc25;
Cs(11,32) = nc19;
Cs(11,33) = nc10;
Cs(11,55) = nc49;
Cs(11,56) = nc40;
Cs(11,57) = nc22;
Cs(11,58) = nc34;
Cs(11,59) = nc16;
Cs(11,60) = nc7;
Cs(11,79) = nc31;
Cs(11,80) = nc13;
Cs(11,81) = nc4;
Cs(11,96) = nc1;
Cs(12,27) = nc58;
Cs(12,34) = nc55;
Cs(12,35) = nc46;
Cs(12,36) = nc28;
Cs(12,40) = nc52;
Cs(12,41) = nc43;
Cs(12,42) = nc25;
Cs(12,46) = nc37;
Cs(12,47) = nc19;
Cs(12,51) = nc10;
Cs(12,61) = nc49;
Cs(12,62) = nc40;
Cs(12,63) = nc22;
Cs(12,67) = nc34;
Cs(12,68) = nc16;
Cs(12,72) = nc7;
Cs(12,82) = nc31;
Cs(12,83) = nc13;
Cs(12,87) = nc4;
Cs(12,97) = nc1;
Cs(13,34) = nc58;
Cs(13,35) = nc55;
Cs(13,36) = nc46;
Cs(13,37) = nc28;
Cs(13,41) = nc52;
Cs(13,42) = nc43;
Cs(13,43) = nc25;
Cs(13,47) = nc37;
Cs(13,48) = nc19;
Cs(13,52) = nc10;
Cs(13,62) = nc49;
Cs(13,63) = nc40;
Cs(13,64) = nc22;
Cs(13,68) = nc34;
Cs(13,69) = nc16;
Cs(13,73) = nc7;
Cs(13,83) = nc31;
Cs(13,84) = nc13;
Cs(13,88) = nc4;
Cs(13,98) = nc1;
Cs(14,35) = nc58;
Cs(14,36) = nc55;
Cs(14,37) = nc46;
Cs(14,38) = nc28;
Cs(14,42) = nc52;
Cs(14,43) = nc43;
Cs(14,44) = nc25;
Cs(14,48) = nc37;
Cs(14,49) = nc19;
Cs(14,53) = nc10;
Cs(14,63) = nc49;
Cs(14,64) = nc40;
Cs(14,65) = nc22;
Cs(14,69) = nc34;
Cs(14,70) = nc16;
Cs(14,74) = nc7;
Cs(14,84) = nc31;
Cs(14,85) = nc13;
Cs(14,89) = nc4;
Cs(14,99) = nc1;
Cs(15,36) = nc58;
Cs(15,37) = nc55;
Cs(15,38) = nc46;
Cs(15,39) = nc28;
Cs(15,43) = nc52;
Cs(15,44) = nc43;
Cs(15,45) = nc25;
Cs(15,49) = nc37;
Cs(15,50) = nc19;
Cs(15,54) = nc10;
Cs(15,64) = nc49;
Cs(15,65) = nc40;
Cs(15,66) = nc22;
Cs(15,70) = nc34;
Cs(15,71) = nc16;
Cs(15,75) = nc7;
Cs(15,85) = nc31;
Cs(15,86) = nc13;
Cs(15,90) = nc4;
Cs(15,100) = nc1;
Cs(16,40) = nc58;
Cs(16,41) = nc55;
Cs(16,42) = nc46;
Cs(16,43) = nc28;
Cs(16,46) = nc52;
Cs(16,47) = nc43;
Cs(16,48) = nc25;
Cs(16,51) = nc37;
Cs(16,52) = nc19;
Cs(16,55) = nc10;
Cs(16,67) = nc49;
Cs(16,68) = nc40;
Cs(16,69) = nc22;
Cs(16,72) = nc34;
Cs(16,73) = nc16;
Cs(16,76) = nc7;
Cs(16,87) = nc31;
Cs(16,88) = nc13;
Cs(16,91) = nc4;
Cs(16,101) = nc1;
Cs(17,41) = nc58;
Cs(17,42) = nc55;
Cs(17,43) = nc46;
Cs(17,44) = nc28;
Cs(17,47) = nc52;
Cs(17,48) = nc43;
Cs(17,49) = nc25;
Cs(17,52) = nc37;
Cs(17,53) = nc19;
Cs(17,56) = nc10;
Cs(17,68) = nc49;
Cs(17,69) = nc40;
Cs(17,70) = nc22;
Cs(17,73) = nc34;
Cs(17,74) = nc16;
Cs(17,77) = nc7;
Cs(17,88) = nc31;
Cs(17,89) = nc13;
Cs(17,92) = nc4;
Cs(17,102) = nc1;
Cs(18,42) = nc58;
Cs(18,43) = nc55;
Cs(18,44) = nc46;
Cs(18,45) = nc28;
Cs(18,48) = nc52;
Cs(18,49) = nc43;
Cs(18,50) = nc25;
Cs(18,53) = nc37;
Cs(18,54) = nc19;
Cs(18,57) = nc10;
Cs(18,69) = nc49;
Cs(18,70) = nc40;
Cs(18,71) = nc22;
Cs(18,74) = nc34;
Cs(18,75) = nc16;
Cs(18,78) = nc7;
Cs(18,89) = nc31;
Cs(18,90) = nc13;
Cs(18,93) = nc4;
Cs(18,103) = nc1;
Cs(19,46) = nc58;
Cs(19,47) = nc55;
Cs(19,48) = nc46;
Cs(19,49) = nc28;
Cs(19,51) = nc52;
Cs(19,52) = nc43;
Cs(19,53) = nc25;
Cs(19,55) = nc37;
Cs(19,56) = nc19;
Cs(19,58) = nc10;
Cs(19,72) = nc49;
Cs(19,73) = nc40;
Cs(19,74) = nc22;
Cs(19,76) = nc34;
Cs(19,77) = nc16;
Cs(19,79) = nc7;
Cs(19,91) = nc31;
Cs(19,92) = nc13;
Cs(19,94) = nc4;
Cs(19,104) = nc1;
Cs(20,47) = nc58;
Cs(20,48) = nc55;
Cs(20,49) = nc46;
Cs(20,50) = nc28;
Cs(20,52) = nc52;
Cs(20,53) = nc43;
Cs(20,54) = nc25;
Cs(20,56) = nc37;
Cs(20,57) = nc19;
Cs(20,59) = nc10;
Cs(20,73) = nc49;
Cs(20,74) = nc40;
Cs(20,75) = nc22;
Cs(20,77) = nc34;
Cs(20,78) = nc16;
Cs(20,80) = nc7;
Cs(20,92) = nc31;
Cs(20,93) = nc13;
Cs(20,95) = nc4;
Cs(20,105) = nc1;
Cs(21,51) = nc58;
Cs(21,52) = nc55;
Cs(21,53) = nc46;
Cs(21,54) = nc28;
Cs(21,55) = nc52;
Cs(21,56) = nc43;
Cs(21,57) = nc25;
Cs(21,58) = nc37;
Cs(21,59) = nc19;
Cs(21,60) = nc10;
Cs(21,76) = nc49;
Cs(21,77) = nc40;
Cs(21,78) = nc22;
Cs(21,79) = nc34;
Cs(21,80) = nc16;
Cs(21,81) = nc7;
Cs(21,94) = nc31;
Cs(21,95) = nc13;
Cs(21,96) = nc4;
Cs(21,106) = nc1;
Cs(22,61) = nc58;
Cs(22,62) = nc55;
Cs(22,63) = nc46;
Cs(22,64) = nc28;
Cs(22,67) = nc52;
Cs(22,68) = nc43;
Cs(22,69) = nc25;
Cs(22,72) = nc37;
Cs(22,73) = nc19;
Cs(22,76) = nc10;
Cs(22,82) = nc49;
Cs(22,83) = nc40;
Cs(22,84) = nc22;
Cs(22,87) = nc34;
Cs(22,88) = nc16;
Cs(22,91) = nc7;
Cs(22,97) = nc31;
Cs(22,98) = nc13;
Cs(22,101) = nc4;
Cs(22,107) = nc1;
Cs(23,62) = nc58;
Cs(23,63) = nc55;
Cs(23,64) = nc46;
Cs(23,65) = nc28;
Cs(23,68) = nc52;
Cs(23,69) = nc43;
Cs(23,70) = nc25;
Cs(23,73) = nc37;
Cs(23,74) = nc19;
Cs(23,77) = nc10;
Cs(23,83) = nc49;
Cs(23,84) = nc40;
Cs(23,85) = nc22;
Cs(23,88) = nc34;
Cs(23,89) = nc16;
Cs(23,92) = nc7;
Cs(23,98) = nc31;
Cs(23,99) = nc13;
Cs(23,102) = nc4;
Cs(23,108) = nc1;
Cs(24,63) = nc58;
Cs(24,64) = nc55;
Cs(24,65) = nc46;
Cs(24,66) = nc28;
Cs(24,69) = nc52;
Cs(24,70) = nc43;
Cs(24,71) = nc25;
Cs(24,74) = nc37;
Cs(24,75) = nc19;
Cs(24,78) = nc10;
Cs(24,84) = nc49;
Cs(24,85) = nc40;
Cs(24,86) = nc22;
Cs(24,89) = nc34;
Cs(24,90) = nc16;
Cs(24,93) = nc7;
Cs(24,99) = nc31;
Cs(24,100) = nc13;
Cs(24,103) = nc4;
Cs(24,109) = nc1;
Cs(25,67) = nc58;
Cs(25,68) = nc55;
Cs(25,69) = nc46;
Cs(25,70) = nc28;
Cs(25,72) = nc52;
Cs(25,73) = nc43;
Cs(25,74) = nc25;
Cs(25,76) = nc37;
Cs(25,77) = nc19;
Cs(25,79) = nc10;
Cs(25,87) = nc49;
Cs(25,88) = nc40;
Cs(25,89) = nc22;
Cs(25,91) = nc34;
Cs(25,92) = nc16;
Cs(25,94) = nc7;
Cs(25,101) = nc31;
Cs(25,102) = nc13;
Cs(25,104) = nc4;
Cs(25,110) = nc1;
Cs(26,68) = nc58;
Cs(26,69) = nc55;
Cs(26,70) = nc46;
Cs(26,71) = nc28;
Cs(26,73) = nc52;
Cs(26,74) = nc43;
Cs(26,75) = nc25;
Cs(26,77) = nc37;
Cs(26,78) = nc19;
Cs(26,80) = nc10;
Cs(26,88) = nc49;
Cs(26,89) = nc40;
Cs(26,90) = nc22;
Cs(26,92) = nc34;
Cs(26,93) = nc16;
Cs(26,95) = nc7;
Cs(26,102) = nc31;
Cs(26,103) = nc13;
Cs(26,105) = nc4;
Cs(26,111) = nc1;
Cs(27,72) = nc58;
Cs(27,73) = nc55;
Cs(27,74) = nc46;
Cs(27,75) = nc28;
Cs(27,76) = nc52;
Cs(27,77) = nc43;
Cs(27,78) = nc25;
Cs(27,79) = nc37;
Cs(27,80) = nc19;
Cs(27,81) = nc10;
Cs(27,91) = nc49;
Cs(27,92) = nc40;
Cs(27,93) = nc22;
Cs(27,94) = nc34;
Cs(27,95) = nc16;
Cs(27,96) = nc7;
Cs(27,104) = nc31;
Cs(27,105) = nc13;
Cs(27,106) = nc4;
Cs(27,112) = nc1;
Cs(28,82) = nc58;
Cs(28,83) = nc55;
Cs(28,84) = nc46;
Cs(28,85) = nc28;
Cs(28,87) = nc52;
Cs(28,88) = nc43;
Cs(28,89) = nc25;
Cs(28,91) = nc37;
Cs(28,92) = nc19;
Cs(28,94) = nc10;
Cs(28,97) = nc49;
Cs(28,98) = nc40;
Cs(28,99) = nc22;
Cs(28,101) = nc34;
Cs(28,102) = nc16;
Cs(28,104) = nc7;
Cs(28,107) = nc31;
Cs(28,108) = nc13;
Cs(28,110) = nc4;
Cs(28,113) = nc1;
Cs(29,87) = nc58;
Cs(29,88) = nc55;
Cs(29,89) = nc46;
Cs(29,90) = nc28;
Cs(29,91) = nc52;
Cs(29,92) = nc43;
Cs(29,93) = nc25;
Cs(29,94) = nc37;
Cs(29,95) = nc19;
Cs(29,96) = nc10;
Cs(29,101) = nc49;
Cs(29,102) = nc40;
Cs(29,103) = nc22;
Cs(29,104) = nc34;
Cs(29,105) = nc16;
Cs(29,106) = nc7;
Cs(29,110) = nc31;
Cs(29,111) = nc13;
Cs(29,112) = nc4;
Cs(29,114) = nc1;
Cs(30,1) = nc59;
Cs(30,2) = nc56;
Cs(30,3) = nc47;
Cs(30,4) = nc29;
Cs(30,7) = nc53;
Cs(30,8) = nc44;
Cs(30,9) = nc26;
Cs(30,12) = nc38;
Cs(30,13) = nc20;
Cs(30,17) = nc11;
Cs(30,27) = nc50;
Cs(30,34) = nc41;
Cs(30,35) = nc23;
Cs(30,40) = nc35;
Cs(30,41) = nc17;
Cs(30,46) = nc8;
Cs(30,61) = nc32;
Cs(30,62) = nc14;
Cs(30,67) = nc5;
Cs(30,82) = nc2;
Cs(31,2) = nc59;
Cs(31,3) = nc56;
Cs(31,4) = nc47;
Cs(31,5) = nc29;
Cs(31,8) = nc53;
Cs(31,9) = nc44;
Cs(31,10) = nc26;
Cs(31,13) = nc38;
Cs(31,14) = nc20;
Cs(31,18) = nc11;
Cs(31,34) = nc50;
Cs(31,35) = nc41;
Cs(31,36) = nc23;
Cs(31,41) = nc35;
Cs(31,42) = nc17;
Cs(31,47) = nc8;
Cs(31,62) = nc32;
Cs(31,63) = nc14;
Cs(31,68) = nc5;
Cs(31,83) = nc2;
Cs(32,3) = nc59;
Cs(32,4) = nc56;
Cs(32,5) = nc47;
Cs(32,6) = nc29;
Cs(32,9) = nc53;
Cs(32,10) = nc44;
Cs(32,11) = nc26;
Cs(32,14) = nc38;
Cs(32,15) = nc20;
Cs(32,19) = nc11;
Cs(32,35) = nc50;
Cs(32,36) = nc41;
Cs(32,37) = nc23;
Cs(32,42) = nc35;
Cs(32,43) = nc17;
Cs(32,48) = nc8;
Cs(32,63) = nc32;
Cs(32,64) = nc14;
Cs(32,69) = nc5;
Cs(32,84) = nc2;
Cs(33,7) = nc59;
Cs(33,8) = nc56;
Cs(33,9) = nc47;
Cs(33,10) = nc29;
Cs(33,12) = nc53;
Cs(33,13) = nc44;
Cs(33,14) = nc26;
Cs(33,17) = nc38;
Cs(33,18) = nc20;
Cs(33,21) = nc11;
Cs(33,40) = nc50;
Cs(33,41) = nc41;
Cs(33,42) = nc23;
Cs(33,46) = nc35;
Cs(33,47) = nc17;
Cs(33,51) = nc8;
Cs(33,67) = nc32;
Cs(33,68) = nc14;
Cs(33,72) = nc5;
Cs(33,87) = nc2;
Cs(34,8) = nc59;
Cs(34,9) = nc56;
Cs(34,10) = nc47;
Cs(34,11) = nc29;
Cs(34,13) = nc53;
Cs(34,14) = nc44;
Cs(34,15) = nc26;
Cs(34,18) = nc38;
Cs(34,19) = nc20;
Cs(34,22) = nc11;
Cs(34,41) = nc50;
Cs(34,42) = nc41;
Cs(34,43) = nc23;
Cs(34,47) = nc35;
Cs(34,48) = nc17;
Cs(34,52) = nc8;
Cs(34,68) = nc32;
Cs(34,69) = nc14;
Cs(34,73) = nc5;
Cs(34,88) = nc2;
Cs(35,12) = nc59;
Cs(35,13) = nc56;
Cs(35,14) = nc47;
Cs(35,15) = nc29;
Cs(35,17) = nc53;
Cs(35,18) = nc44;
Cs(35,19) = nc26;
Cs(35,21) = nc38;
Cs(35,22) = nc20;
Cs(35,24) = nc11;
Cs(35,46) = nc50;
Cs(35,47) = nc41;
Cs(35,48) = nc23;
Cs(35,51) = nc35;
Cs(35,52) = nc17;
Cs(35,55) = nc8;
Cs(35,72) = nc32;
Cs(35,73) = nc14;
Cs(35,76) = nc5;
Cs(35,91) = nc2;
Cs(36,13) = nc59;
Cs(36,14) = nc56;
Cs(36,15) = nc47;
Cs(36,16) = nc29;
Cs(36,18) = nc53;
Cs(36,19) = nc44;
Cs(36,20) = nc26;
Cs(36,22) = nc38;
Cs(36,23) = nc20;
Cs(36,25) = nc11;
Cs(36,47) = nc50;
Cs(36,48) = nc41;
Cs(36,49) = nc23;
Cs(36,52) = nc35;
Cs(36,53) = nc17;
Cs(36,56) = nc8;
Cs(36,73) = nc32;
Cs(36,74) = nc14;
Cs(36,77) = nc5;
Cs(36,92) = nc2;
Cs(37,14) = nc59;
Cs(37,15) = nc56;
Cs(37,16) = nc47;
Cs(37,19) = nc53;
Cs(37,20) = nc44;
Cs(37,23) = nc38;
Cs(37,28) = nc29;
Cs(37,29) = nc26;
Cs(37,30) = nc20;
Cs(37,31) = nc11;
Cs(37,48) = nc50;
Cs(37,49) = nc41;
Cs(37,50) = nc23;
Cs(37,53) = nc35;
Cs(37,54) = nc17;
Cs(37,57) = nc8;
Cs(37,74) = nc32;
Cs(37,75) = nc14;
Cs(37,78) = nc5;
Cs(37,93) = nc2;
Cs(38,17) = nc59;
Cs(38,18) = nc56;
Cs(38,19) = nc47;
Cs(38,20) = nc29;
Cs(38,21) = nc53;
Cs(38,22) = nc44;
Cs(38,23) = nc26;
Cs(38,24) = nc38;
Cs(38,25) = nc20;
Cs(38,26) = nc11;
Cs(38,51) = nc50;
Cs(38,52) = nc41;
Cs(38,53) = nc23;
Cs(38,55) = nc35;
Cs(38,56) = nc17;
Cs(38,58) = nc8;
Cs(38,76) = nc32;
Cs(38,77) = nc14;
Cs(38,79) = nc5;
Cs(38,94) = nc2;
Cs(39,18) = nc59;
Cs(39,19) = nc56;
Cs(39,20) = nc47;
Cs(39,22) = nc53;
Cs(39,23) = nc44;
Cs(39,25) = nc38;
Cs(39,29) = nc29;
Cs(39,30) = nc26;
Cs(39,31) = nc20;
Cs(39,32) = nc11;
Cs(39,52) = nc50;
Cs(39,53) = nc41;
Cs(39,54) = nc23;
Cs(39,56) = nc35;
Cs(39,57) = nc17;
Cs(39,59) = nc8;
Cs(39,77) = nc32;
Cs(39,78) = nc14;
Cs(39,80) = nc5;
Cs(39,95) = nc2;
Cs(40,21) = nc59;
Cs(40,22) = nc56;
Cs(40,23) = nc47;
Cs(40,24) = nc53;
Cs(40,25) = nc44;
Cs(40,26) = nc38;
Cs(40,30) = nc29;
Cs(40,31) = nc26;
Cs(40,32) = nc20;
Cs(40,33) = nc11;
Cs(40,55) = nc50;
Cs(40,56) = nc41;
Cs(40,57) = nc23;
Cs(40,58) = nc35;
Cs(40,59) = nc17;
Cs(40,60) = nc8;
Cs(40,79) = nc32;
Cs(40,80) = nc14;
Cs(40,81) = nc5;
Cs(40,96) = nc2;
Cs(41,27) = nc59;
Cs(41,34) = nc56;
Cs(41,35) = nc47;
Cs(41,36) = nc29;
Cs(41,40) = nc53;
Cs(41,41) = nc44;
Cs(41,42) = nc26;
Cs(41,46) = nc38;
Cs(41,47) = nc20;
Cs(41,51) = nc11;
Cs(41,61) = nc50;
Cs(41,62) = nc41;
Cs(41,63) = nc23;
Cs(41,67) = nc35;
Cs(41,68) = nc17;
Cs(41,72) = nc8;
Cs(41,82) = nc32;
Cs(41,83) = nc14;
Cs(41,87) = nc5;
Cs(41,97) = nc2;
Cs(42,34) = nc59;
Cs(42,35) = nc56;
Cs(42,36) = nc47;
Cs(42,37) = nc29;
Cs(42,41) = nc53;
Cs(42,42) = nc44;
Cs(42,43) = nc26;
Cs(42,47) = nc38;
Cs(42,48) = nc20;
Cs(42,52) = nc11;
Cs(42,62) = nc50;
Cs(42,63) = nc41;
Cs(42,64) = nc23;
Cs(42,68) = nc35;
Cs(42,69) = nc17;
Cs(42,73) = nc8;
Cs(42,83) = nc32;
Cs(42,84) = nc14;
Cs(42,88) = nc5;
Cs(42,98) = nc2;
Cs(43,35) = nc59;
Cs(43,36) = nc56;
Cs(43,37) = nc47;
Cs(43,38) = nc29;
Cs(43,42) = nc53;
Cs(43,43) = nc44;
Cs(43,44) = nc26;
Cs(43,48) = nc38;
Cs(43,49) = nc20;
Cs(43,53) = nc11;
Cs(43,63) = nc50;
Cs(43,64) = nc41;
Cs(43,65) = nc23;
Cs(43,69) = nc35;
Cs(43,70) = nc17;
Cs(43,74) = nc8;
Cs(43,84) = nc32;
Cs(43,85) = nc14;
Cs(43,89) = nc5;
Cs(43,99) = nc2;
Cs(44,36) = nc59;
Cs(44,37) = nc56;
Cs(44,38) = nc47;
Cs(44,39) = nc29;
Cs(44,43) = nc53;
Cs(44,44) = nc44;
Cs(44,45) = nc26;
Cs(44,49) = nc38;
Cs(44,50) = nc20;
Cs(44,54) = nc11;
Cs(44,64) = nc50;
Cs(44,65) = nc41;
Cs(44,66) = nc23;
Cs(44,70) = nc35;
Cs(44,71) = nc17;
Cs(44,75) = nc8;
Cs(44,85) = nc32;
Cs(44,86) = nc14;
Cs(44,90) = nc5;
Cs(44,100) = nc2;
Cs(45,40) = nc59;
Cs(45,41) = nc56;
Cs(45,42) = nc47;
Cs(45,43) = nc29;
Cs(45,46) = nc53;
Cs(45,47) = nc44;
Cs(45,48) = nc26;
Cs(45,51) = nc38;
Cs(45,52) = nc20;
Cs(45,55) = nc11;
Cs(45,67) = nc50;
Cs(45,68) = nc41;
Cs(45,69) = nc23;
Cs(45,72) = nc35;
Cs(45,73) = nc17;
Cs(45,76) = nc8;
Cs(45,87) = nc32;
Cs(45,88) = nc14;
Cs(45,91) = nc5;
Cs(45,101) = nc2;
Cs(46,41) = nc59;
Cs(46,42) = nc56;
Cs(46,43) = nc47;
Cs(46,44) = nc29;
Cs(46,47) = nc53;
Cs(46,48) = nc44;
Cs(46,49) = nc26;
Cs(46,52) = nc38;
Cs(46,53) = nc20;
Cs(46,56) = nc11;
Cs(46,68) = nc50;
Cs(46,69) = nc41;
Cs(46,70) = nc23;
Cs(46,73) = nc35;
Cs(46,74) = nc17;
Cs(46,77) = nc8;
Cs(46,88) = nc32;
Cs(46,89) = nc14;
Cs(46,92) = nc5;
Cs(46,102) = nc2;
Cs(47,42) = nc59;
Cs(47,43) = nc56;
Cs(47,44) = nc47;
Cs(47,45) = nc29;
Cs(47,48) = nc53;
Cs(47,49) = nc44;
Cs(47,50) = nc26;
Cs(47,53) = nc38;
Cs(47,54) = nc20;
Cs(47,57) = nc11;
Cs(47,69) = nc50;
Cs(47,70) = nc41;
Cs(47,71) = nc23;
Cs(47,74) = nc35;
Cs(47,75) = nc17;
Cs(47,78) = nc8;
Cs(47,89) = nc32;
Cs(47,90) = nc14;
Cs(47,93) = nc5;
Cs(47,103) = nc2;
Cs(48,46) = nc59;
Cs(48,47) = nc56;
Cs(48,48) = nc47;
Cs(48,49) = nc29;
Cs(48,51) = nc53;
Cs(48,52) = nc44;
Cs(48,53) = nc26;
Cs(48,55) = nc38;
Cs(48,56) = nc20;
Cs(48,58) = nc11;
Cs(48,72) = nc50;
Cs(48,73) = nc41;
Cs(48,74) = nc23;
Cs(48,76) = nc35;
Cs(48,77) = nc17;
Cs(48,79) = nc8;
Cs(48,91) = nc32;
Cs(48,92) = nc14;
Cs(48,94) = nc5;
Cs(48,104) = nc2;
Cs(49,47) = nc59;
Cs(49,48) = nc56;
Cs(49,49) = nc47;
Cs(49,50) = nc29;
Cs(49,52) = nc53;
Cs(49,53) = nc44;
Cs(49,54) = nc26;
Cs(49,56) = nc38;
Cs(49,57) = nc20;
Cs(49,59) = nc11;
Cs(49,73) = nc50;
Cs(49,74) = nc41;
Cs(49,75) = nc23;
Cs(49,77) = nc35;
Cs(49,78) = nc17;
Cs(49,80) = nc8;
Cs(49,92) = nc32;
Cs(49,93) = nc14;
Cs(49,95) = nc5;
Cs(49,105) = nc2;
Cs(50,51) = nc59;
Cs(50,52) = nc56;
Cs(50,53) = nc47;
Cs(50,54) = nc29;
Cs(50,55) = nc53;
Cs(50,56) = nc44;
Cs(50,57) = nc26;
Cs(50,58) = nc38;
Cs(50,59) = nc20;
Cs(50,60) = nc11;
Cs(50,76) = nc50;
Cs(50,77) = nc41;
Cs(50,78) = nc23;
Cs(50,79) = nc35;
Cs(50,80) = nc17;
Cs(50,81) = nc8;
Cs(50,94) = nc32;
Cs(50,95) = nc14;
Cs(50,96) = nc5;
Cs(50,106) = nc2;
Cs(51,61) = nc59;
Cs(51,62) = nc56;
Cs(51,63) = nc47;
Cs(51,64) = nc29;
Cs(51,67) = nc53;
Cs(51,68) = nc44;
Cs(51,69) = nc26;
Cs(51,72) = nc38;
Cs(51,73) = nc20;
Cs(51,76) = nc11;
Cs(51,82) = nc50;
Cs(51,83) = nc41;
Cs(51,84) = nc23;
Cs(51,87) = nc35;
Cs(51,88) = nc17;
Cs(51,91) = nc8;
Cs(51,97) = nc32;
Cs(51,98) = nc14;
Cs(51,101) = nc5;
Cs(51,107) = nc2;
Cs(52,62) = nc59;
Cs(52,63) = nc56;
Cs(52,64) = nc47;
Cs(52,65) = nc29;
Cs(52,68) = nc53;
Cs(52,69) = nc44;
Cs(52,70) = nc26;
Cs(52,73) = nc38;
Cs(52,74) = nc20;
Cs(52,77) = nc11;
Cs(52,83) = nc50;
Cs(52,84) = nc41;
Cs(52,85) = nc23;
Cs(52,88) = nc35;
Cs(52,89) = nc17;
Cs(52,92) = nc8;
Cs(52,98) = nc32;
Cs(52,99) = nc14;
Cs(52,102) = nc5;
Cs(52,108) = nc2;
Cs(53,63) = nc59;
Cs(53,64) = nc56;
Cs(53,65) = nc47;
Cs(53,66) = nc29;
Cs(53,69) = nc53;
Cs(53,70) = nc44;
Cs(53,71) = nc26;
Cs(53,74) = nc38;
Cs(53,75) = nc20;
Cs(53,78) = nc11;
Cs(53,84) = nc50;
Cs(53,85) = nc41;
Cs(53,86) = nc23;
Cs(53,89) = nc35;
Cs(53,90) = nc17;
Cs(53,93) = nc8;
Cs(53,99) = nc32;
Cs(53,100) = nc14;
Cs(53,103) = nc5;
Cs(53,109) = nc2;
Cs(54,67) = nc59;
Cs(54,68) = nc56;
Cs(54,69) = nc47;
Cs(54,70) = nc29;
Cs(54,72) = nc53;
Cs(54,73) = nc44;
Cs(54,74) = nc26;
Cs(54,76) = nc38;
Cs(54,77) = nc20;
Cs(54,79) = nc11;
Cs(54,87) = nc50;
Cs(54,88) = nc41;
Cs(54,89) = nc23;
Cs(54,91) = nc35;
Cs(54,92) = nc17;
Cs(54,94) = nc8;
Cs(54,101) = nc32;
Cs(54,102) = nc14;
Cs(54,104) = nc5;
Cs(54,110) = nc2;
Cs(55,68) = nc59;
Cs(55,69) = nc56;
Cs(55,70) = nc47;
Cs(55,71) = nc29;
Cs(55,73) = nc53;
Cs(55,74) = nc44;
Cs(55,75) = nc26;
Cs(55,77) = nc38;
Cs(55,78) = nc20;
Cs(55,80) = nc11;
Cs(55,88) = nc50;
Cs(55,89) = nc41;
Cs(55,90) = nc23;
Cs(55,92) = nc35;
Cs(55,93) = nc17;
Cs(55,95) = nc8;
Cs(55,102) = nc32;
Cs(55,103) = nc14;
Cs(55,105) = nc5;
Cs(55,111) = nc2;
Cs(56,72) = nc59;
Cs(56,73) = nc56;
Cs(56,74) = nc47;
Cs(56,75) = nc29;
Cs(56,76) = nc53;
Cs(56,77) = nc44;
Cs(56,78) = nc26;
Cs(56,79) = nc38;
Cs(56,80) = nc20;
Cs(56,81) = nc11;
Cs(56,91) = nc50;
Cs(56,92) = nc41;
Cs(56,93) = nc23;
Cs(56,94) = nc35;
Cs(56,95) = nc17;
Cs(56,96) = nc8;
Cs(56,104) = nc32;
Cs(56,105) = nc14;
Cs(56,106) = nc5;
Cs(56,112) = nc2;
Cs(57,82) = nc59;
Cs(57,83) = nc56;
Cs(57,84) = nc47;
Cs(57,85) = nc29;
Cs(57,87) = nc53;
Cs(57,88) = nc44;
Cs(57,89) = nc26;
Cs(57,91) = nc38;
Cs(57,92) = nc20;
Cs(57,94) = nc11;
Cs(57,97) = nc50;
Cs(57,98) = nc41;
Cs(57,99) = nc23;
Cs(57,101) = nc35;
Cs(57,102) = nc17;
Cs(57,104) = nc8;
Cs(57,107) = nc32;
Cs(57,108) = nc14;
Cs(57,110) = nc5;
Cs(57,113) = nc2;
Cs(58,87) = nc59;
Cs(58,88) = nc56;
Cs(58,89) = nc47;
Cs(58,90) = nc29;
Cs(58,91) = nc53;
Cs(58,92) = nc44;
Cs(58,93) = nc26;
Cs(58,94) = nc38;
Cs(58,95) = nc20;
Cs(58,96) = nc11;
Cs(58,101) = nc50;
Cs(58,102) = nc41;
Cs(58,103) = nc23;
Cs(58,104) = nc35;
Cs(58,105) = nc17;
Cs(58,106) = nc8;
Cs(58,110) = nc32;
Cs(58,111) = nc14;
Cs(58,112) = nc5;
Cs(58,114) = nc2;
Cs(59,1) = nc60;
Cs(59,2) = nc57;
Cs(59,3) = nc48;
Cs(59,4) = nc30;
Cs(59,7) = nc54;
Cs(59,8) = nc45;
Cs(59,9) = nc27;
Cs(59,12) = nc39;
Cs(59,13) = nc21;
Cs(59,17) = nc12;
Cs(59,27) = nc51;
Cs(59,34) = nc42;
Cs(59,35) = nc24;
Cs(59,40) = nc36;
Cs(59,41) = nc18;
Cs(59,46) = nc9;
Cs(59,61) = nc33;
Cs(59,62) = nc15;
Cs(59,67) = nc6;
Cs(59,82) = nc3;
Cs(60,2) = nc60;
Cs(60,3) = nc57;
Cs(60,4) = nc48;
Cs(60,5) = nc30;
Cs(60,8) = nc54;
Cs(60,9) = nc45;
Cs(60,10) = nc27;
Cs(60,13) = nc39;
Cs(60,14) = nc21;
Cs(60,18) = nc12;
Cs(60,34) = nc51;
Cs(60,35) = nc42;
Cs(60,36) = nc24;
Cs(60,41) = nc36;
Cs(60,42) = nc18;
Cs(60,47) = nc9;
Cs(60,62) = nc33;
Cs(60,63) = nc15;
Cs(60,68) = nc6;
Cs(60,83) = nc3;
Cs(61,3) = nc60;
Cs(61,4) = nc57;
Cs(61,5) = nc48;
Cs(61,6) = nc30;
Cs(61,9) = nc54;
Cs(61,10) = nc45;
Cs(61,11) = nc27;
Cs(61,14) = nc39;
Cs(61,15) = nc21;
Cs(61,19) = nc12;
Cs(61,35) = nc51;
Cs(61,36) = nc42;
Cs(61,37) = nc24;
Cs(61,42) = nc36;
Cs(61,43) = nc18;
Cs(61,48) = nc9;
Cs(61,63) = nc33;
Cs(61,64) = nc15;
Cs(61,69) = nc6;
Cs(61,84) = nc3;
Cs(62,7) = nc60;
Cs(62,8) = nc57;
Cs(62,9) = nc48;
Cs(62,10) = nc30;
Cs(62,12) = nc54;
Cs(62,13) = nc45;
Cs(62,14) = nc27;
Cs(62,17) = nc39;
Cs(62,18) = nc21;
Cs(62,21) = nc12;
Cs(62,40) = nc51;
Cs(62,41) = nc42;
Cs(62,42) = nc24;
Cs(62,46) = nc36;
Cs(62,47) = nc18;
Cs(62,51) = nc9;
Cs(62,67) = nc33;
Cs(62,68) = nc15;
Cs(62,72) = nc6;
Cs(62,87) = nc3;
Cs(63,8) = nc60;
Cs(63,9) = nc57;
Cs(63,10) = nc48;
Cs(63,11) = nc30;
Cs(63,13) = nc54;
Cs(63,14) = nc45;
Cs(63,15) = nc27;
Cs(63,18) = nc39;
Cs(63,19) = nc21;
Cs(63,22) = nc12;
Cs(63,41) = nc51;
Cs(63,42) = nc42;
Cs(63,43) = nc24;
Cs(63,47) = nc36;
Cs(63,48) = nc18;
Cs(63,52) = nc9;
Cs(63,68) = nc33;
Cs(63,69) = nc15;
Cs(63,73) = nc6;
Cs(63,88) = nc3;
Cs(64,12) = nc60;
Cs(64,13) = nc57;
Cs(64,14) = nc48;
Cs(64,15) = nc30;
Cs(64,17) = nc54;
Cs(64,18) = nc45;
Cs(64,19) = nc27;
Cs(64,21) = nc39;
Cs(64,22) = nc21;
Cs(64,24) = nc12;
Cs(64,46) = nc51;
Cs(64,47) = nc42;
Cs(64,48) = nc24;
Cs(64,51) = nc36;
Cs(64,52) = nc18;
Cs(64,55) = nc9;
Cs(64,72) = nc33;
Cs(64,73) = nc15;
Cs(64,76) = nc6;
Cs(64,91) = nc3;
Cs(65,13) = nc60;
Cs(65,14) = nc57;
Cs(65,15) = nc48;
Cs(65,16) = nc30;
Cs(65,18) = nc54;
Cs(65,19) = nc45;
Cs(65,20) = nc27;
Cs(65,22) = nc39;
Cs(65,23) = nc21;
Cs(65,25) = nc12;
Cs(65,47) = nc51;
Cs(65,48) = nc42;
Cs(65,49) = nc24;
Cs(65,52) = nc36;
Cs(65,53) = nc18;
Cs(65,56) = nc9;
Cs(65,73) = nc33;
Cs(65,74) = nc15;
Cs(65,77) = nc6;
Cs(65,92) = nc3;
Cs(66,14) = nc60;
Cs(66,15) = nc57;
Cs(66,16) = nc48;
Cs(66,19) = nc54;
Cs(66,20) = nc45;
Cs(66,23) = nc39;
Cs(66,28) = nc30;
Cs(66,29) = nc27;
Cs(66,30) = nc21;
Cs(66,31) = nc12;
Cs(66,48) = nc51;
Cs(66,49) = nc42;
Cs(66,50) = nc24;
Cs(66,53) = nc36;
Cs(66,54) = nc18;
Cs(66,57) = nc9;
Cs(66,74) = nc33;
Cs(66,75) = nc15;
Cs(66,78) = nc6;
Cs(66,93) = nc3;
Cs(67,17) = nc60;
Cs(67,18) = nc57;
Cs(67,19) = nc48;
Cs(67,20) = nc30;
Cs(67,21) = nc54;
Cs(67,22) = nc45;
Cs(67,23) = nc27;
Cs(67,24) = nc39;
Cs(67,25) = nc21;
Cs(67,26) = nc12;
Cs(67,51) = nc51;
Cs(67,52) = nc42;
Cs(67,53) = nc24;
Cs(67,55) = nc36;
Cs(67,56) = nc18;
Cs(67,58) = nc9;
Cs(67,76) = nc33;
Cs(67,77) = nc15;
Cs(67,79) = nc6;
Cs(67,94) = nc3;
Cs(68,18) = nc60;
Cs(68,19) = nc57;
Cs(68,20) = nc48;
Cs(68,22) = nc54;
Cs(68,23) = nc45;
Cs(68,25) = nc39;
Cs(68,29) = nc30;
Cs(68,30) = nc27;
Cs(68,31) = nc21;
Cs(68,32) = nc12;
Cs(68,52) = nc51;
Cs(68,53) = nc42;
Cs(68,54) = nc24;
Cs(68,56) = nc36;
Cs(68,57) = nc18;
Cs(68,59) = nc9;
Cs(68,77) = nc33;
Cs(68,78) = nc15;
Cs(68,80) = nc6;
Cs(68,95) = nc3;
Cs(69,21) = nc60;
Cs(69,22) = nc57;
Cs(69,23) = nc48;
Cs(69,24) = nc54;
Cs(69,25) = nc45;
Cs(69,26) = nc39;
Cs(69,30) = nc30;
Cs(69,31) = nc27;
Cs(69,32) = nc21;
Cs(69,33) = nc12;
Cs(69,55) = nc51;
Cs(69,56) = nc42;
Cs(69,57) = nc24;
Cs(69,58) = nc36;
Cs(69,59) = nc18;
Cs(69,60) = nc9;
Cs(69,79) = nc33;
Cs(69,80) = nc15;
Cs(69,81) = nc6;
Cs(69,96) = nc3;
Cs(70,27) = nc60;
Cs(70,34) = nc57;
Cs(70,35) = nc48;
Cs(70,36) = nc30;
Cs(70,40) = nc54;
Cs(70,41) = nc45;
Cs(70,42) = nc27;
Cs(70,46) = nc39;
Cs(70,47) = nc21;
Cs(70,51) = nc12;
Cs(70,61) = nc51;
Cs(70,62) = nc42;
Cs(70,63) = nc24;
Cs(70,67) = nc36;
Cs(70,68) = nc18;
Cs(70,72) = nc9;
Cs(70,82) = nc33;
Cs(70,83) = nc15;
Cs(70,87) = nc6;
Cs(70,97) = nc3;
Cs(71,34) = nc60;
Cs(71,35) = nc57;
Cs(71,36) = nc48;
Cs(71,37) = nc30;
Cs(71,41) = nc54;
Cs(71,42) = nc45;
Cs(71,43) = nc27;
Cs(71,47) = nc39;
Cs(71,48) = nc21;
Cs(71,52) = nc12;
Cs(71,62) = nc51;
Cs(71,63) = nc42;
Cs(71,64) = nc24;
Cs(71,68) = nc36;
Cs(71,69) = nc18;
Cs(71,73) = nc9;
Cs(71,83) = nc33;
Cs(71,84) = nc15;
Cs(71,88) = nc6;
Cs(71,98) = nc3;
Cs(72,35) = nc60;
Cs(72,36) = nc57;
Cs(72,37) = nc48;
Cs(72,38) = nc30;
Cs(72,42) = nc54;
Cs(72,43) = nc45;
Cs(72,44) = nc27;
Cs(72,48) = nc39;
Cs(72,49) = nc21;
Cs(72,53) = nc12;
Cs(72,63) = nc51;
Cs(72,64) = nc42;
Cs(72,65) = nc24;
Cs(72,69) = nc36;
Cs(72,70) = nc18;
Cs(72,74) = nc9;
Cs(72,84) = nc33;
Cs(72,85) = nc15;
Cs(72,89) = nc6;
Cs(72,99) = nc3;
Cs(73,36) = nc60;
Cs(73,37) = nc57;
Cs(73,38) = nc48;
Cs(73,39) = nc30;
Cs(73,43) = nc54;
Cs(73,44) = nc45;
Cs(73,45) = nc27;
Cs(73,49) = nc39;
Cs(73,50) = nc21;
Cs(73,54) = nc12;
Cs(73,64) = nc51;
Cs(73,65) = nc42;
Cs(73,66) = nc24;
Cs(73,70) = nc36;
Cs(73,71) = nc18;
Cs(73,75) = nc9;
Cs(73,85) = nc33;
Cs(73,86) = nc15;
Cs(73,90) = nc6;
Cs(73,100) = nc3;
Cs(74,40) = nc60;
Cs(74,41) = nc57;
Cs(74,42) = nc48;
Cs(74,43) = nc30;
Cs(74,46) = nc54;
Cs(74,47) = nc45;
Cs(74,48) = nc27;
Cs(74,51) = nc39;
Cs(74,52) = nc21;
Cs(74,55) = nc12;
Cs(74,67) = nc51;
Cs(74,68) = nc42;
Cs(74,69) = nc24;
Cs(74,72) = nc36;
Cs(74,73) = nc18;
Cs(74,76) = nc9;
Cs(74,87) = nc33;
Cs(74,88) = nc15;
Cs(74,91) = nc6;
Cs(74,101) = nc3;
Cs(75,41) = nc60;
Cs(75,42) = nc57;
Cs(75,43) = nc48;
Cs(75,44) = nc30;
Cs(75,47) = nc54;
Cs(75,48) = nc45;
Cs(75,49) = nc27;
Cs(75,52) = nc39;
Cs(75,53) = nc21;
Cs(75,56) = nc12;
Cs(75,68) = nc51;
Cs(75,69) = nc42;
Cs(75,70) = nc24;
Cs(75,73) = nc36;
Cs(75,74) = nc18;
Cs(75,77) = nc9;
Cs(75,88) = nc33;
Cs(75,89) = nc15;
Cs(75,92) = nc6;
Cs(75,102) = nc3;
Cs(76,42) = nc60;
Cs(76,43) = nc57;
Cs(76,44) = nc48;
Cs(76,45) = nc30;
Cs(76,48) = nc54;
Cs(76,49) = nc45;
Cs(76,50) = nc27;
Cs(76,53) = nc39;
Cs(76,54) = nc21;
Cs(76,57) = nc12;
Cs(76,69) = nc51;
Cs(76,70) = nc42;
Cs(76,71) = nc24;
Cs(76,74) = nc36;
Cs(76,75) = nc18;
Cs(76,78) = nc9;
Cs(76,89) = nc33;
Cs(76,90) = nc15;
Cs(76,93) = nc6;
Cs(76,103) = nc3;
Cs(77,46) = nc60;
Cs(77,47) = nc57;
Cs(77,48) = nc48;
Cs(77,49) = nc30;
Cs(77,51) = nc54;
Cs(77,52) = nc45;
Cs(77,53) = nc27;
Cs(77,55) = nc39;
Cs(77,56) = nc21;
Cs(77,58) = nc12;
Cs(77,72) = nc51;
Cs(77,73) = nc42;
Cs(77,74) = nc24;
Cs(77,76) = nc36;
Cs(77,77) = nc18;
Cs(77,79) = nc9;
Cs(77,91) = nc33;
Cs(77,92) = nc15;
Cs(77,94) = nc6;
Cs(77,104) = nc3;
Cs(78,47) = nc60;
Cs(78,48) = nc57;
Cs(78,49) = nc48;
Cs(78,50) = nc30;
Cs(78,52) = nc54;
Cs(78,53) = nc45;
Cs(78,54) = nc27;
Cs(78,56) = nc39;
Cs(78,57) = nc21;
Cs(78,59) = nc12;
Cs(78,73) = nc51;
Cs(78,74) = nc42;
Cs(78,75) = nc24;
Cs(78,77) = nc36;
Cs(78,78) = nc18;
Cs(78,80) = nc9;
Cs(78,92) = nc33;
Cs(78,93) = nc15;
Cs(78,95) = nc6;
Cs(78,105) = nc3;
Cs(79,51) = nc60;
Cs(79,52) = nc57;
Cs(79,53) = nc48;
Cs(79,54) = nc30;
Cs(79,55) = nc54;
Cs(79,56) = nc45;
Cs(79,57) = nc27;
Cs(79,58) = nc39;
Cs(79,59) = nc21;
Cs(79,60) = nc12;
Cs(79,76) = nc51;
Cs(79,77) = nc42;
Cs(79,78) = nc24;
Cs(79,79) = nc36;
Cs(79,80) = nc18;
Cs(79,81) = nc9;
Cs(79,94) = nc33;
Cs(79,95) = nc15;
Cs(79,96) = nc6;
Cs(79,106) = nc3;
Cs(80,61) = nc60;
Cs(80,62) = nc57;
Cs(80,63) = nc48;
Cs(80,64) = nc30;
Cs(80,67) = nc54;
Cs(80,68) = nc45;
Cs(80,69) = nc27;
Cs(80,72) = nc39;
Cs(80,73) = nc21;
Cs(80,76) = nc12;
Cs(80,82) = nc51;
Cs(80,83) = nc42;
Cs(80,84) = nc24;
Cs(80,87) = nc36;
Cs(80,88) = nc18;
Cs(80,91) = nc9;
Cs(80,97) = nc33;
Cs(80,98) = nc15;
Cs(80,101) = nc6;
Cs(80,107) = nc3;
Cs(81,62) = nc60;
Cs(81,63) = nc57;
Cs(81,64) = nc48;
Cs(81,65) = nc30;
Cs(81,68) = nc54;
Cs(81,69) = nc45;
Cs(81,70) = nc27;
Cs(81,73) = nc39;
Cs(81,74) = nc21;
Cs(81,77) = nc12;
Cs(81,83) = nc51;
Cs(81,84) = nc42;
Cs(81,85) = nc24;
Cs(81,88) = nc36;
Cs(81,89) = nc18;
Cs(81,92) = nc9;
Cs(81,98) = nc33;
Cs(81,99) = nc15;
Cs(81,102) = nc6;
Cs(81,108) = nc3;
Cs(82,63) = nc60;
Cs(82,64) = nc57;
Cs(82,65) = nc48;
Cs(82,66) = nc30;
Cs(82,69) = nc54;
Cs(82,70) = nc45;
Cs(82,71) = nc27;
Cs(82,74) = nc39;
Cs(82,75) = nc21;
Cs(82,78) = nc12;
Cs(82,84) = nc51;
Cs(82,85) = nc42;
Cs(82,86) = nc24;
Cs(82,89) = nc36;
Cs(82,90) = nc18;
Cs(82,93) = nc9;
Cs(82,99) = nc33;
Cs(82,100) = nc15;
Cs(82,103) = nc6;
Cs(82,109) = nc3;
Cs(83,67) = nc60;
Cs(83,68) = nc57;
Cs(83,69) = nc48;
Cs(83,70) = nc30;
Cs(83,72) = nc54;
Cs(83,73) = nc45;
Cs(83,74) = nc27;
Cs(83,76) = nc39;
Cs(83,77) = nc21;
Cs(83,79) = nc12;
Cs(83,87) = nc51;
Cs(83,88) = nc42;
Cs(83,89) = nc24;
Cs(83,91) = nc36;
Cs(83,92) = nc18;
Cs(83,94) = nc9;
Cs(83,101) = nc33;
Cs(83,102) = nc15;
Cs(83,104) = nc6;
Cs(83,110) = nc3;
Cs(84,68) = nc60;
Cs(84,69) = nc57;
Cs(84,70) = nc48;
Cs(84,71) = nc30;
Cs(84,73) = nc54;
Cs(84,74) = nc45;
Cs(84,75) = nc27;
Cs(84,77) = nc39;
Cs(84,78) = nc21;
Cs(84,80) = nc12;
Cs(84,88) = nc51;
Cs(84,89) = nc42;
Cs(84,90) = nc24;
Cs(84,92) = nc36;
Cs(84,93) = nc18;
Cs(84,95) = nc9;
Cs(84,102) = nc33;
Cs(84,103) = nc15;
Cs(84,105) = nc6;
Cs(84,111) = nc3;
Cs(85,72) = nc60;
Cs(85,73) = nc57;
Cs(85,74) = nc48;
Cs(85,75) = nc30;
Cs(85,76) = nc54;
Cs(85,77) = nc45;
Cs(85,78) = nc27;
Cs(85,79) = nc39;
Cs(85,80) = nc21;
Cs(85,81) = nc12;
Cs(85,91) = nc51;
Cs(85,92) = nc42;
Cs(85,93) = nc24;
Cs(85,94) = nc36;
Cs(85,95) = nc18;
Cs(85,96) = nc9;
Cs(85,104) = nc33;
Cs(85,105) = nc15;
Cs(85,106) = nc6;
Cs(85,112) = nc3;
Cs(86,82) = nc60;
Cs(86,83) = nc57;
Cs(86,84) = nc48;
Cs(86,85) = nc30;
Cs(86,87) = nc54;
Cs(86,88) = nc45;
Cs(86,89) = nc27;
Cs(86,91) = nc39;
Cs(86,92) = nc21;
Cs(86,94) = nc12;
Cs(86,97) = nc51;
Cs(86,98) = nc42;
Cs(86,99) = nc24;
Cs(86,101) = nc36;
Cs(86,102) = nc18;
Cs(86,104) = nc9;
Cs(86,107) = nc33;
Cs(86,108) = nc15;
Cs(86,110) = nc6;
Cs(86,113) = nc3;
Cs(87,87) = nc60;
Cs(87,88) = nc57;
Cs(87,89) = nc48;
Cs(87,90) = nc30;
Cs(87,91) = nc54;
Cs(87,92) = nc45;
Cs(87,93) = nc27;
Cs(87,94) = nc39;
Cs(87,95) = nc21;
Cs(87,96) = nc12;
Cs(87,101) = nc51;
Cs(87,102) = nc42;
Cs(87,103) = nc24;
Cs(87,104) = nc36;
Cs(87,105) = nc18;
Cs(87,106) = nc9;
Cs(87,110) = nc33;
Cs(87,111) = nc15;
Cs(87,112) = nc6;
Cs(87,114) = nc3;
Cs(88,27) = 1;
Cs(88,115) = -1;
Cs(89,34) = 1;
Cs(89,116) = -1;
Cs(90,35) = 1;
Cs(90,117) = -1;
Cs(91,36) = 1;
Cs(91,118) = -1;
Cs(92,37) = 1;
Cs(92,119) = -1;
Cs(93,38) = 1;
Cs(93,120) = -1;
Cs(94,40) = 1;
Cs(94,121) = -1;
Cs(95,41) = 1;
Cs(95,122) = -1;
Cs(96,42) = 1;
Cs(96,123) = -1;
Cs(97,43) = 1;
Cs(97,124) = -1;
Cs(98,44) = 1;
Cs(98,125) = -1;
Cs(99,46) = 1;
Cs(99,126) = -1;
Cs(100,47) = 1;
Cs(100,127) = -1;
Cs(101,48) = 1;
Cs(101,128) = -1;
Cs(102,49) = 1;
Cs(102,129) = -1;
Cs(103,50) = 1;
Cs(103,130) = -1;
Cs(104,51) = 1;
Cs(104,131) = -1;
Cs(105,52) = 1;
Cs(105,132) = -1;
Cs(106,53) = 1;
Cs(106,133) = -1;
Cs(107,54) = 1;
Cs(107,134) = -1;
Cs(108,55) = 1;
Cs(108,135) = -1;
Cs(109,56) = 1;
Cs(109,136) = -1;
Cs(110,57) = 1;
Cs(110,137) = -1;
Cs(111,58) = 1;
Cs(111,138) = -1;
Cs(112,59) = 1;
Cs(112,139) = -1;
Cs(113,60) = 1;
Cs(113,140) = -1;
Cs(114,61) = 1;
Cs(114,141) = -1;
solForm = [[-1] [-1] [-1]; [0] [0] [1]; [0] [0] [0]; [0] [0] [0]; [0] [0] [0]; [0] [0] [0]; [0] [1] [0]; [0] [0] [0]; [0] [0] [0]; [0] [0] [0]; [0] [0] [0]; [0] [0] [0]; [0] [0] [0]; [0] [0] [0]; [0] [0] [0]; [0] [0] [0]; [0] [0] [0]; [0] [0] [0]; [0] [0] [0]; [0] [0] [0]; [0] [0] [0]; [0] [0] [0]; [0] [0] [0]; [0] [0] [0]; [0] [0] [0]; [0] [0] [0]; [1] [0] [0];];

allCss = mat2cell(Cs, 114, ones(1,2) * 114);
noOfVars = 3; 
C0 = allCss{1}; 
C1 = allCss{2};
A1 = C0(1:end-27,1:27);
A2 = C0(1:end-27,28:end);
B1 = C0(end-27+1 : end,1:27);
B2 = C0(end-27+1 : end, 28:end);
X = B1 - B2 * (A2 \ A1); 
[V,D] = eig(X);
EValues = diag(D); 
EVectors = V;
good = ~(isinf(EValues) | isnan(EValues));
EValues = EValues(good);
EVectors = EVectors(:,good);

% The eigen values and eigen vectors have been now extracted
PEPsolutions=[];
nonInfEValuesInd = ~isinf(EValues);
NinfEValues = EValues(nonInfEValuesInd);
NinfEVectors = EVectors(:,nonInfEValuesInd);
noOfEvalues = length(NinfEValues);
sizeOfEvectors = length(solForm);

% We basically, then iterate through all of the received eigenvalues and
% then try to remove those that do not satisfy the criterion for the
% corresponding eigenvectors to have a form that is the same as that of the
% monomial vector.
% In fact we also remove those eigenvalues and eigenvectors which give us
% solutions that have infinity value for atleast one variable.
g = 1;
% monstocheck(rowstorem1)=[];
% monstocheck(rowstorem2)=[];
allvarsextracted = sum(abs(solForm));
if length(find(allvarsextracted==0)) == 0
    
    for i = 1:noOfEvalues
        otherVarValues = NinfEVectors(:,i);
        for k = 1:noOfVars
            sols(i,k) = 1;
            for j = 1:sizeOfEvectors
                sols(i,k) = sols(i,k) * otherVarValues(j) ^ solForm(j,k);
            end
        end
        %     sols(i,1) = NinfEValues(i);
        sols(i,k + 1) = NinfEValues(i);
        PEPsolutions(g,:) = sols(i,:);
        g = g+1;
        
    end

else
    PEPsolutions = [];
end
end
