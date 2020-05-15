function sols_data = solver_problem_l2_3view_triang(data)
C = setup_elimination_template(data);
C0 = C(:,1:250);
C1 = C(:,251:end);
C1 = C0 \ C1;
RR = [-C1(end-28:end,:);eye(31)];
AM_ind = [39,1,4,5,2,6,3,7,8,9,10,11,12,13,14,15,16,17,40,18,19,20,21,22,23,24,25,26,27,28,29];
AM = RR(AM_ind,:);
[V,D] = eig(AM);
V = V ./ (ones(size(V,1),1)*V(1,:));
sols(1,:) = V(2,:);
sols(2,:) = V(5,:);
sols(3,:) = diag(D).';
sols(4,:) = V(12,:);
sols(5,:) = V(16,:);
sols(6,:) = V(19,:);
sols(7,:) = V(23,:);
sols(8,:) = V(29,:);
sols_data = {C,sols}; 

% Action =  x3
% Quotient ring basis (V) = 1,x1,x1*x4,x1*x8,x2,x2*x4,x2*x6,x2*x8,x2*x8^2,x3,x3*x6,x4,x4^2,x4*x6,x4*x6^2,x5,x5*x6,x5*x8,x6,x6^2,x6*x8,x6*x8^2,x7,x7^2,x7^3,x7^2*x8,x7*x8,x7*x8^2,x8,x8^2,x8^3,
% Available monomials (RR*V) = x1*x3,x2*x3,x2*x3*x6,x1*x3*x4,x1*x3*x8,x2*x3*x4,x2*x3*x8,x2*x3*x8^2,x3^2,x3^2*x6,x3*x4,x3*x4^2,x3*x4*x6,x3*x4*x6^2,x3*x5,x3*x5*x6,x3*x5*x8,x3*x6^2,x3*x6*x8,x3*x6*x8^2,x3*x7,x3*x7^2,x3*x7^3,x3*x7^2*x8,x3*x7*x8,x3*x7*x8^2,x3*x8,x3*x8^2,x3*x8^3,1,x1,x1*x4,x1*x8,x2,x2*x4,x2*x6,x2*x8,x2*x8^2,x3,x3*x6,x4,x4^2,x4*x6,x4*x6^2,x5,x5*x6,x5*x8,x6,x6^2,x6*x8,x6*x8^2,x7,x7^2,x7^3,x7^2*x8,x7*x8,x7*x8^2,x8,x8^2,x8^3,
function [coeffs,coeffs_ind] = compute_coeffs(data)
coeffs(1) = data(1);
coeffs(2) = data(4);
coeffs(3) = 2;
coeffs(4) = data(7);
coeffs(5) = -2*data(19);
coeffs(6) = data(2);
coeffs(7) = data(5);
coeffs(8) = data(8);
coeffs(9) = -2*data(20);
coeffs(10) = data(10);
coeffs(11) = data(13);
coeffs(12) = data(3);
coeffs(13) = data(16);
coeffs(14) = -2*data(21);
coeffs(15) = data(11);
coeffs(16) = data(14);
coeffs(17) = data(6);
coeffs(18) = data(17);
coeffs(19) = -2*data(22);
coeffs(20) = data(12);
coeffs(21) = -2*data(23);
coeffs(22) = data(15);
coeffs(23) = -2*data(24);
coeffs(24) = data(9);
coeffs(25) = data(18);
coeffs_ind = [1,10,6,10,2,1,2,15,7,15,6,7,1,1,11,1,2,2,16,16,6,7,7,1,1,2,2,10,11,15,10,16,11,15,10,11,16,15,16,10,10,1,10,6,1,6,10,2,1,2,...
15,7,2,6,7,7,15,6,1,2,1,11,1,6,6,7,6,11,7,1,2,2,16,2,6,7,7,7,16,10,6,1,15,7,2,11,6,1,16,7,2,11,10,10,16,15,15,10,15,11,...
10,11,10,11,10,16,15,16,15,16,15,11,16,10,11,11,11,16,15,16,16,10,10,11,11,3,10,3,3,10,15,3,3,15,3,3,11,3,3,16,1,2,6,1,2,7,6,7,1,1,...
2,1,7,6,6,1,6,1,1,2,6,1,2,2,7,6,7,7,6,2,2,1,1,7,6,1,7,7,2,6,6,7,6,2,1,7,1,2,7,6,2,7,7,6,1,2,7,2,1,10,...
6,10,2,1,2,15,7,15,6,7,6,1,10,7,2,15,1,1,11,13,6,1,1,2,2,16,18,18,16,6,7,7,7,6,1,11,2,7,2,16,4,1,2,20,8,20,6,7,12,10,...
4,1,2,22,22,8,6,7,1,1,3,1,2,2,1,6,2,3,2,7,6,1,10,7,2,15,6,1,11,7,2,16,4,2,1,8,7,6,17,17,2,7,10,11,10,11,10,11,1,15,...
10,11,16,6,15,10,16,11,15,16,2,15,16,7,15,10,11,16,15,16,3,10,15,3,10,15,1,10,10,6,10,15,10,10,2,15,15,3,15,11,10,10,15,3,10,16,15,3,3,11,...
16,3,11,16,1,11,20,3,22,16,10,11,11,20,3,3,22,16,15,10,11,11,11,16,16,20,22,3,11,15,16,16,10,10,3,10,15,3,3,20,11,16,10,15,3,22,11,16,3,2,...
1,6,1,13,13,7,1,6,6,6,1,6,2,7,1,2,18,18,2,7,7,6,7,2,7,1,6,2,7,4,3,1,2,20,8,4,8,3,20,6,7,12,10,6,1,10,11,17,12,17,...
15,7,2,15,16,8,3,4,2,1,22,4,8,3,8,22,7,6,12,17,11,12,10,1,6,11,11,12,17,16,17,15,2,7,16,16,18,15,16,2,7,10,11,1,15,16,2,15,7,16,...
13,15,10,10,11,10,13,18,18,15,15,16,10,15,22,3,20,10,15,13,18,16,15,18,11,16,16,22,20,22,11,16,15,10,20,20,3,22,16,11,22,2,16,7,15,18,15,16,3,3,...
1,3,3,6,3,1,3,3,6,3,3,3,2,3,3,7,3,3,3,3,2,3,3,7,1,10,3,6,10,3,10,3,3,2,15,3,3,7,15,3,3,10,15,3,3,15,3,1,11,3,...
3,2,16,3,3,7,16,21,23,18,3,3,16,14,19,20,3,3,3,15,2,1,12,17,4,7,2,1,6,12,17,8,7,6,3,4,1,6,1,2,13,12,8,13,12,8,4,6,1,6,...
7,6,17,8,3,12,4,1,6,12,8,3,17,4,7,2,1,2,18,12,17,8,4,18,3,17,8,2,7,6,7,7,8,17,17,12,12,4,13,6,1,2,7,12,8,4,3,17,18,7,...
2,4,3,2,1,20,8,20,3,7,6,12,6,10,1,17,12,17,8,4,20,7,15,2,3,4,2,22,25,1,25,22,3,8,7,6,12,8,13,12,4,11,1,6,12,17,8,4,18,17,...
22,16,2,7,24,12,17,20,3,3,24,22,12,17,3,3,4,1,2,8,6,7,3,19,17,4,8,3,17,18,7,2,8,4,3,20,8,4,22,3,13,18,3,3,13,18,3,18,20,3,...
22,10,11,13,13,20,22,18,3,4,16,15,2,18,13,18,20,22,8,16,15,7,18,3,20,18,3,13,13,12,22,15,11,3,16,18,10,13,20,22,3,18,17,16,15,4,3,3,15,10,...
20,5,9,8,3,3,15,10,20,14,19,13,18,17,22,3,20,10,15,15,16,23,3,15,10,9,4,3,3,16,11,22,5,9,8,16,11,3,3,22,14,19,13,18,17,20,22,3,22,18,...
11,16,16,16,15,13,18,21,23,3,3,20,16,11,15,10,21,25,13,18,23,22,16,3,11,3,8,4,1,2,25,25,4,8,8,6,7,3,18,17,4,12,17,7,18,15,8,16,2,24,...
9,5,12,17,8,4,10,15,20,22,3,9,5,24,17,12,20,4,8,22,11,16,22,3,20,22,3,4,18,13,17,15,18,16,21,25,13,18,20,22,3,10,15,20,14,19,4,3,3,14,...
19,8,3,4,3,3,20,8,20,3,3,14,19,9,5,23,17,20,3,15,3,4,3,3,22,17,8,4,12,2,1,25,3,25,24,17,12,4,8,3,7,6,8,19,5,9,24,25,17,8,...
4,18,3,12,7,17,2,24,9,5,17,12,8,20,4,23,3,25,9,5,24,21,22,17,4,23,8,3,12,19,5,25,24,22,14,13,3,17,18,15,9,18,16,20,4,3,3,1,13,14,...
9,21,5,12,8,23,13,3,6,3,3,3,6,11,3,2,3,1,3,11,10,1,13,3,3,7,6,13,20,3,22,11,10,6,13,10,6,11,19,9,21,5,12,23,3,13,21,23,13,3,...
3,11,9,21,19,17,14,5,12,23,13,3,3,18,3,3,21,23,13,18,3,3,11,16,3,3,11,16,12,3,3,10,20,20,3,3,10,12,3,10,15,20,10,11,14,22,3,3,11,14,...
19,13,12,20,3,13,22,11,11,16,22,11,10,13,18,16,11,10,13,11,11,14,5,24,12,12,13,3,4,17,6,8,1,14,12,4,3,12,13,6,1,12,12,1,6,13,10,11,1,6,...
13,17,4,12,8,12,6,13,10,11,1,1,11,6,10,19,5,14,9,13,20,22,3,18,25,12,13,10,11,18,13,12,10,13,11,13,10,11,25,19,9,21,24,14,5,23,19,9,23,14,...
21,4,25,3,14,9,21,23,19,4,3,3,2,18,21,23,18,9,22,13,20,4,25,3,24,25,19,14,5,23,21,3,8,9,14,19,9,5,21,23,17,8,18,3,3,7,8,22,3,3,...
5,21,18,13,9,22,23,3,20,8,25,20,8,22,3,24,14,21,23,5,19,13,12,9,25,3,25,14,9,21,5,12,23,3,11,22,13,3,24,25,9,5,23,14,18,3,17,19,21,14,...
9,5,21,23,17,18,3,25,14,19,9,5,21,23,22,17,3,3,16,18,14,19,22,3,3,16,24,19,14,9,5,20,3,23,14,19,20,3,3,24,14,19,9,5,21,23,18,13,3,15,...
20,10,22,24,19,5,23,22,9,25,21,14,3,14,19,3,22,24,19,14,5,21,25,18,9,22,23,16,22,11,20,13,3,23,18,13,20,25,22,11,16,22,3,24,25,14,19,9,17,8,...
21,12,23,4,5,19,5,21,14,24,25,12,17,8,23,4,9,24,17,4,12,8,18,20,13,22,4,8,13,18,5,24,12,17,8,25,23,20,22,9,4,21,4,22,8,20,14,19,9,21,...
24,18,25,20,13,5,22,23,5,21,24,13,18,9,20,25,23,22,25,20,22];
function C = setup_elimination_template(data)
[coeffs,coeffs_ind] = compute_coeffs(data);
C_ind = [1,235,251,252,501,503,504,735,751,752,753,754,1005,1014,1235,1256,1257,1264,1485,1502,1506,1507,1514,1758,1774,2009,2024,2260,2285,2510,2511,2535,2536,2761,2762,2763,2786,3012,3013,3255,3288,3515,3739,3765,3766,3768,3862,4015,4017,4028,...
4239,4265,4266,4267,4268,4278,4362,4519,4529,4531,4611,4739,4770,4772,4779,4781,4861,4862,5019,5021,5030,5111,5239,5270,5271,5272,5280,5361,5362,5508,5523,5575,5759,5773,5825,6008,6025,6026,6259,6275,6276,6527,6596,6767,6777,6778,6846,7029,7031,7096,...
7103,7108,7136,7267,7271,7278,7280,7346,7353,7358,7386,7529,7531,7532,7604,7636,7771,7780,7782,7854,7886,8033,8134,8284,8384,8510,8539,8753,8761,8787,8789,9004,9012,9037,9255,9285,9289,9507,9513,9537,9790,9845,10040,10041,10042,10095,10291,10292,10543,10598,...
10794,10849,11044,11045,11099,11140,11296,11392,11547,11600,11797,11798,11799,11848,11850,12048,12049,12098,12300,12349,12351,12394,12395,12545,12551,12554,12599,12601,12640,12644,12645,12796,12802,12892,12897,13050,13053,13145,13301,13303,13304,13395,13552,13555,13556,13647,13805,13806,14057,14234,...
14307,14308,14557,14559,14560,14734,14807,14808,14809,14810,15061,15062,15153,15311,15312,15403,15563,15609,15734,15735,15814,15906,16065,16066,16109,16234,16235,16252,16308,16315,16316,16359,16564,16567,16568,16653,16656,16817,16818,16903,17001,17069,17070,17235,17251,17252,17319,17320,17501,17571,...
17764,17822,17823,17985,18002,18014,18072,18073,18324,18416,18508,18584,18826,18916,19077,19078,19084,19259,19327,19328,19579,19580,19676,19829,19830,19926,20081,20082,20176,20331,20332,20426,20524,20583,20637,20774,20833,20887,21009,21024,21085,21244,21336,21400,21587,21651,21838,21841,21902,22086,...
22089,22092,22150,22152,22337,22340,22401,22404,22588,22591,22652,22839,22842,22902,23090,23093,23094,23154,23343,23344,23516,23540,23595,23768,23791,23792,24097,24227,24293,24347,24475,24544,24711,24797,24847,24850,24977,25017,25101,25102,25144,25146,25211,25278,25303,25352,25396,25519,25520,25540,...
25595,25772,25791,25792,26105,26227,26260,26279,26285,26294,26356,26410,26461,26511,26521,26531,26536,26601,26606,26607,26610,26644,26646,26660,26711,26762,26763,26780,26803,26857,26860,26896,27063,27159,27277,27319,27320,27353,27532,27538,27569,27570,27572,27573,27608,27788,27822,27823,27854,28113,...
28117,28186,28200,28239,28362,28363,28364,28367,28437,28450,28615,28618,28866,28936,28942,28950,28989,29112,29114,29116,29187,29192,29200,29365,29368,29369,29476,29619,29726,29765,29825,29870,29871,29989,30015,30016,30018,30023,30112,30120,30121,30265,30324,30372,30373,30455,30469,30515,30517,30528,...
30576,30622,30623,30705,30719,30769,30776,30861,30874,30875,30989,31020,31022,31025,31111,31112,31124,31125,31279,31281,31324,31361,31376,31377,31429,31430,31455,31521,31530,31576,31611,31626,31627,31679,31680,31705,31759,31878,31934,31935,31995,32129,32130,32193,32379,32380,32443,32631,32693,32694,...
32846,32863,32867,32882,32883,32959,33017,33028,33096,33116,33132,33133,33192,33209,33277,33283,33346,33370,33371,33521,33530,33616,33635,33636,33692,33701,33709,33846,33853,33858,33870,33871,33874,33875,33886,34032,34034,34104,34124,34125,34136,34388,34456,34457,34497,34634,34639,34736,34790,34890,...
34891,35041,35045,35141,35392,35393,35543,35546,35643,35800,35804,35845,35891,36042,36051,36141,36294,36394,36556,36600,36648,36799,36805,36898,37149,37224,37336,37399,37408,37587,37729,37812,37838,37899,37974,38061,38089,38149,38158,38309,38340,38405,38479,38560,38593,38655,38906,38907,38974,39068,...
39091,39157,39224,39317,39342,39407,39408,39512,39513,39537,39566,39594,39655,39753,39754,39787,39820,39852,39857,39987,40161,40197,40290,40345,40348,40411,40412,40413,40447,40541,40542,40598,40662,40663,40834,40849,40914,40917,40918,40960,40984,41043,41045,41058,41098,41099,41140,41164,41165,41168,...
41210,41232,41294,41296,41324,41349,41392,41448,41473,41547,41550,41577,41600,41645,41667,41668,41695,41733,41734,41798,41799,41801,41804,41808,41828,41848,41895,41915,41918,41945,41982,41983,42052,42099,42101,42144,42145,42147,42153,42169,42170,42198,42223,42303,42305,42306,42326,42395,42403,42419,...
42420,42557,42580,42671,42672,42734,42807,42808,42829,42921,42922,43057,43173,43174,43178,43307,43309,43310,43311,43312,43403,43423,43424,43428,43582,43609,43675,43734,43735,43741,43752,43808,43831,43859,43925,43991,44063,44064,44071,44109,44156,44174,44177,44228,44315,44316,44317,44318,44321,44359,...
44403,44424,44427,44478,44501,44569,44570,44571,44626,44719,44764,44821,44822,44823,44930,45085,45166,45181,45182,45416,45431,45432,45494,45509,45576,45577,45578,45583,45666,45676,45683,45743,45829,45830,45878,45926,46081,46082,46176,46184,46290,46345,46364,46436,46541,46542,46687,46794,46837,46867,...
46901,46938,46946,46961,47047,47088,47091,47100,47119,47152,47189,47212,47215,47227,47298,47299,47339,47342,47402,47440,47463,47465,47475,47476,47590,47601,47613,47644,47646,47652,47654,47688,47691,47692,47696,47711,47750,47803,47843,47844,47866,47896,47902,47941,48000,48097,48123,48129,48161,48197,...
48227,48266,48268,48347,48372,48381,48412,48413,48475,48517,48528,48559,48560,48597,48602,48621,48646,48695,48717,48733,48738,48777,48882,48921,48922,49019,49105,49127,49130,49161,49197,49227,49270,49272,49355,49412,49413,49429,49444,49475,49521,49530,49565,49566,49605,49607,49610,49624,49646,49659,...
49695,49699,49717,49733,49740,49819,49820,49853,49858,49883,49885,49909,49921,49922,49925,49991,50032,50038,50072,50073,50104,50159,50175,50201,50241,50435,50436,50450,50452,50453,50489,50612,50614,50687,50700,50702,50703,50745,50826,50866,50869,50942,50950,50954,50955,50970,50976,50981,50996,51015,...
51023,51075,51120,51121,51122,51123,51181,51182,51205,51219,51247,51275,51276,51361,51374,51375,51376,51377,51429,51430,51431,51432,51455,51456,51629,51630,51638,51693,51866,51942,51943,51958,51959,51999,52033,52096,52120,52121,52132,52133,52139,52202,52203,52209,52290,52345,52391,52415,52417,52541,...
52542,52641,52732,52899,52928,52962,52974,53149,53158,53173,53213,53309,53310,53311,53312,53352,53399,53405,53421,53464,53500,53657,53677,53689,53724,53911,53917,53918,53947,53966,53980,53984,53993,54058,54098,54162,54163,54165,54168,54183,54216,54230,54232,54326,54327,54328,54395,54403,54418,54419,...
54420,54424,54432,54445,54468,54483,54492,54557,54579,54580,54671,54672,54673,54674,54678,54719,54720,54821,54831,54832,54859,54876,54924,54925,54927,54930,54978,54981,54991,55116,55119,55146,55152,55191,55192,55195,55203,55215,55217,55221,55226,55233,55248,55250,55393,55447,55448,55472,55474,55543,...
55545,55586,55640,55641,55643,55650,55658,55663,55722,55723,55814,55900,55907,55908,56047,56143,56147,56148,56365,56439,56462,56465,56477,56548,56552,56643,56648,56793,56836,56868,56900,56940,56963,56965,56975,57131,57193,57194,57294,57296,57337,57392,57393,57401,57414,57479,57510,57535,57539,57563,...
57651,57729,57802,57840,57851,57893,57894,57897,57898,57904,57905,57945,57960,57979,58053,58101,58261,58286,58287,58289,58315,58404,58405,58479,58506,58536,58537,58539,58899,58922,58938,58964,58979,59039,59069,59106,59237,59347,59370,59414,59460,59461,59467,59488,59505,59539,59572,59660,59737,59779,...
59781,59813,59855,59856,59875,59909,59910,59914,59949,59960,59961,59967,59990,60029,60031,60113,60117,60135,60136,60201,60209,60324,60334,60349,60414,60418,60424,60431,60448,60460,60468,60473,60492,60508,60574,60584,60637,60666,60676,60683,60743,60758,60774,60835,60994,61008,61128,61184,61185,61245,...
61324,61363,61365,61367,61368,61450,61454,61455,61470,61481,61496,61638,61706,61707,61747,61863,61865,61867,61868,61914,61938,61946,61952,61960,61961,61965,61967,61971,61998,62113,62117,62193,62208,62209,62249,62384,62389,62486,62714,62716,62718,62721,62722,62730,62742,62748,62911,62917,62939,62947,...
62962,62972,62974,62992,63047,63050,63088,63091,63100,63148,63161,63170,63222,63224,63379,63380,63411,63436,63439,63447,63462,63465,63477,63496,63641,63658,63662,63663,63665,63690,63713,63718,63722,63732,63798,63799,63801,63804,63839,63842,63891,63898,63908,63912,63919,63972,64157,64158,64190,64228,...
64364,64381,64412,64413,64437,64440,64444,64454,64463,64465,64475,64631,64693,64694,64707,64893,64914,64938,64946,64948,64960,64964,64972,64973,64979,64980,65039,65063,65064,65106,65156,65157,65160,65196,65214,65229,65237,65241,65398,65405,65419,65420,65441,65445,65464,65466,65472,65483,65500,65553,...
65555,65556,65593,65594,65648,65655,65733,65787,65815,65816,65817,65818,65857,65860,65905,65907,65925,65941,65964,65987,66006,66007,66037,66073,66110,66237,66399,66421,66422,66423,66428,66464,66471,66488,66569,66570,66737,66738,66740,66847,66870,66871,66872,66873,66882,66883,66921,66922,66958,66966,...
66967,66980,66988,67157,67175,67177,67199,67214,67228,67237,67240,67241,67248,67322,67323,67449,67487,67605,67624,67625,67627,67635,67659,67675,67679,67699,67701,67716,67717,67730,67740,67741,67749,67784,67874,67875,67885,67886,67951,67952,67953,67959,67986,68168,68174,68181,68182,68183,68216,68218,...
68220,68230,68231,68242,68243,68333,68335,68378,68387,68416,68426,68431,68432,68433,68434,68493,68494,68524,68583,68585,68637,68744,68833,68878,68887,68934,68935,68995,69181,69182,69185,69200,69202,69203,69204,69205,69206,69220,69231,69245,69246,69247,69388,69456,69457,69497,69702,69703,69704,69708,...
69715,69716,69717,69721,69730,69746,69748,69749,69888,69889,69943,69952,69953,69957,69958,69959,69986,69999,70134,70139,70236];
C = zeros(250,281);
C(C_ind) = coeffs(coeffs_ind);

