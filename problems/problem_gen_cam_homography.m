function[eqsHandler, cfg] = problem_gen_cam_homography()
%% Configuring the solver
cfg = retrieve_solver_cfg();
% If the polynomial system is huge, this function can be
% extracted out in a separate file and the file can be loaded here to
% obtain the input polynomial system.s
eqsHandler = @retrieve_eqs;
end
%% The polynomial system of the solver
% The parameters are the variables and the coefficients. The variables have
% to be labelled as 'a1', 'a2', ... and the coefficients are labelled as
% 'c1', 'c2', ... 
% function eqs = retrieve_eqs(f1, nz, ny, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14, d15, d16, d17, d18, d19, d20, d21, d22, d23, d24, d25, d26, d27, d28, d29, d30, d31, d32, d33, d34, d35, d36, d37, d38, d39, d40, d41, d42, d43, d44, d45, d46, d47, d48, d49, d50, d51, d52, d53, d54, d55, d56, d57, d58, d59, d60, d61, d62, d63, d64, d65, d66, d67, d68, d69, d70, d71)
function eqs = retrieve_eqs(a34,a44,t2,t3,h33, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17, c18)
eqs = [10772*a34^2*a44^2*t2^2+10772*a34^2*a44^2*t3^2+15299*a34^2*a44^2*t2+26996*a34^2*a44^2*t3+3338*a34^2*a44*h33*t2+4979*a34^2*a44*h33*t3+12593*a34^2*a44*t2^2+12593*a34^2*a44*t3^2+11734*a34^2*a44^2+10193*a34^2*a44*h33+11645*a34^2*a44*t2+25887*a34^2*a44*t3+1618*a34^2*h33^2+21067*a34^2*h33*t2+17705*a34^2*h33*t3+17716*a34^2*t2^2+17716*a34^2*t3^2+17850*a34*a44*t2^2+17850*a34*a44*t3^2+29683*a34^2*a44+21143*a34^2*h33+5173*a34^2*t2+18116*a34^2*t3+5154*a34*a44*t2+11698*a34*a44*t3+6953*a34*h33*t2+20626*a34*h33*t3+21185*a34*t2^2+21185*a34*t3^2+21273*a34^2+2655*a34*a44+13012*a34*h33+29591*a34*t2+11750*a34*t3+15170*t2^2+15170*t3^2+24081*a34+8086*t2+26427*t3+29825, 20715*a34^2*a44^2*t2^2+20715*a34^2*a44^2*t3^2+22236*a34^2*a44^2*t2+6843*a34^2*a44^2*t3+452*a34^2*a44*h33*t2+25396*a34^2*a44*h33*t3+17755*a34^2*a44*t2^2+17755*a34^2*a44*t3^2+14602*a34^2*a44^2+8790*a34^2*a44*h33+20977*a34^2*a44*t2+26206*a34^2*a44*t3+14664*a34^2*h33^2+24483*a34^2*h33*t2+17557*a34^2*h33*t3+10440*a34^2*t2^2+10440*a34^2*t3^2+19969*a34*a44*t2^2+19969*a34*a44*t3^2+21925*a34^2*a44+2684*a34^2*h33+4178*a34^2*t2+6059*a34^2*t3+27692*a34*a44*t2+23101*a34*a44*t3+20906*a34*h33*t2+6211*a34*h33*t3+6837*a34*t2^2+6837*a34*t3^2+12940*a34^2+10695*a34*a44+23969*a34*h33+7824*a34*t2+3792*a34*t3+1391*t2^2+1391*t3^2+2936*a34+565*t2+5553*t3+20468, 26234*a34^2*a44^2*t2^2+26234*a34^2*a44^2*t3^2+13642*a34^2*a44^2*t2+6351*a34^2*a44^2*t3+5684*a34^2*a44*h33*t2+6116*a34^2*a44*h33*t3+29934*a34^2*a44*t2^2+29934*a34^2*a44*t3^2+25000*a34^2*a44^2+28810*a34^2*a44*h33+6428*a34^2*a44*t2+3134*a34^2*a44*t3+29229*a34^2*h33^2+29639*a34^2*h33*t2+22447*a34^2*h33*t3+6698*a34^2*t2^2+6698*a34^2*t3^2+20588*a34*a44*t2^2+20588*a34*a44*t3^2+12734*a34^2*a44+5603*a34^2*h33+578*a34^2*t2+1587*a34^2*t3+27182*a34*a44*t2+22916*a34*a44*t3+4584*a34*h33*t2+17401*a34*h33*t3+581*a34*t2^2+581*a34*t3^2+24839*a34^2+1589*a34*a44+18961*a34*h33+1520*a34*t2+10831*a34*t3+744*t2^2+744*t3^2+18720*a34+19354*t2+29407*t3+24609, 27256*a34^2*a44^2*t2^2+27256*a34^2*a44^2*t3^2+5619*a34^2*a44^2*t2+13363*a34^2*a44^2*t3+25710*a34^2*a44*h33*t2+12424*a34^2*a44*h33*t3+10363*a34^2*a44*t2^2+10363*a34^2*a44*t3^2+216*a34^2*a44^2+20267*a34^2*a44*h33+10681*a34^2*a44*t2+9536*a34^2*a44*t3+13504*a34^2*h33^2+20041*a34^2*h33*t2+3232*a34^2*h33*t3+19410*a34^2*t2^2+19410*a34^2*t3^2+24110*a34*a44*t2^2+24110*a34*a44*t3^2+1765*a34^2*a44+27351*a34^2*h33+9253*a34^2*t2+22144*a34^2*t3+8107*a34*a44*t2+3641*a34*a44*t3+13869*a34*h33*t2+8472*a34*h33*t3+22493*a34*t2^2+22493*a34*t3^2+29497*a34^2+26287*a34*a44+17164*a34*h33+14556*a34*t2+11494*a34*t3+16108*t2^2+16108*t3^2+12544*a34+7004*t2+15547*t3+944, 29385*a34^2*a44^2*t3+9277*a34^2*a44*h33*t3+18043*a34^2*a44^2+6624*a34^2*a44*h33+4276*a34^2*a44*t3+22*a34^2*h33^2+10546*a34^2*h33*t3+28790*a34^2*a44+3091*a34^2*h33+23847*a34^2*t3+18170*a34*a44*t3+2463*a34*h33*t3+20968*a34^2+11643*a34*a44+13995*a34*h33+23342*a34*t3+21828*a34+941*t3+28372, 1477*a34^2*a44^2*t2+3455*a34^2*a44*h33*t2+11688*a34^2*a44^2+7749*a34^2*a44*h33+29850*a34^2*a44*t2+21282*a34^2*h33^2+5261*a34^2*h33*t2+10557*a34^2*a44+2267*a34^2*h33+7470*a34^2*t2+22586*a34*a44*t2+90*a34*h33*t2+15254*a34^2+16173*a34*a44+15044*a34*h33+25891*a34*t2+13470*a34+17704*t2+7848, 4578*a34^2*a44^2*t2^2+27256*a34^2*a44^2*t3^2+10200*a34^2*a44^2*t2+13363*a34^2*a44^2*t3+801*a34^2*a44*h33*t2+12424*a34^2*a44*h33*t3+28030*a34^2*a44*t2^2+10363*a34^2*a44*t3^2+14202*a34^2*a44^2+14127*a34^2*a44*h33+3080*a34^2*a44*t2+9536*a34^2*a44*t3+15827*a34^2*h33^2+8471*a34^2*h33*t2+3232*a34^2*h33*t3+29013*a34^2*t2^2+19410*a34^2*t3^2+17539*a34*a44*t2^2+24110*a34*a44*t3^2+14942*a34^2*a44+3181*a34^2*h33+13735*a34^2*t2+22144*a34^2*t3+8540*a34*a44*t2+3641*a34*a44*t3+23269*a34*h33*t2+8472*a34*h33*t3+96*a34*t2^2+22493*a34*t3^2+26971*a34^2+10212*a34*a44+9439*a34*h33+19452*a34*t2+11494*a34*t3+27900*t2^2+16108*t3^2+7642*a34+3619*t2+15547*t3+29772, 8998*a34^2*a44^2*t2^2+8998*a34^2*a44^2*t3^2+381*a34^2*a44^2*t2+229*a34^2*a44^2*t3+948*a34^2*a44*h33*t2+6606*a34^2*a44*h33*t3+28223*a34^2*a44*t2^2+28223*a34^2*a44*t3^2+11464*a34^2*a44^2+16165*a34^2*a44*h33+6001*a34^2*a44*t2+6903*a34^2*a44*t3+8814*a34^2*h33^2+13897*a34^2*h33*t2+7000*a34^2*h33*t3+20611*a34^2*t2^2+20611*a34^2*t3^2+26022*a34*a44*t2^2+26022*a34*a44*t3^2+6033*a34^2*a44+8683*a34^2*h33+15085*a34^2*t2+10151*a34^2*t3+28047*a34*a44*t2+10071*a34*a44*t3+23573*a34*h33*t2+6928*a34*h33*t3+24920*a34*t2^2+24920*a34*t3^2+23902*a34^2+5993*a34*a44+24543*a34*h33+22186*a34*t2+27105*a34*t3+20615*t2^2+20615*t3^2+6167*a34+13550*t2+4500*t3+24066];
end

%% The polynomial system
function cfg = retrieve_solver_cfg()
cfg = {};

% Number of coefficients, labelled as c1, c2, c3,...
cfg.numOfCoeff = 18;

% Number of variables, labeled as a1, a2, a3,...
cfg.numOfVars = 5;

% The index i of the selected variable, ai in the extra polynomial,
% ai - lambda.
% If set to -1, all variables will be tested one
% by one
cfg.hiddenVarNum = 1;

% (1) Either the size of polynomial combinations to be tested.
cfg.sizeOfCombs = [3];
% (2) Or the specific polynomial combination to be tested.
cfg.polyComb=[];
% (3) Or if both are given, the polycomb takes precedence over sizeofcombs.

% The number of rows to be GJ eliminated to obtained a reduced input
% polynomial system as an input to the generator.
cfg.noOfRowsToReduce = 0;

% The heuristic size of the template. There is no theoretical backing, to the best of our
% knowledge, governing the smallest template that can be generated.
% Hence one can start with a larger size and try to test by reducing the size
% of the template.
cfg.heurisiticTemplatesize = 650;
end