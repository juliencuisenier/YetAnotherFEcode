// Gmsh project created on Tue Jun 21 15:23:40 2022
SetFactory("OpenCASCADE");
Mesh.MshFileVersion = 2.2;
Mesh.SubdivisionAlgorithm = 2;

h=2;


//CENTERS
Point(1) = {0,0,0,h};
Point(2) = {0,0,1,h};
Point(3) = {0,0,0.1,h};
Point(4) = {0,0,0.9,h};


//TOP BLADE

//Bottom part
Point(5) = {0.15, 0.476969, 0, h};
Point(6) = {0.15, 0.476969, 1, h};
Point(7) = {-0.15, 0.476969, 0, h};
Point(8) = {-0.15, 0.476969, 1, h};

Point(9) = {0.15, 0.476969, 0.1, h};
Point(10) = {0.15, 0.476969, 0.9, h};
Point(11) = {-0.15, 0.476969, 0.1, h};
Point(12) = {-0.15, 0.476969, 0.9, h};

//Top part

Point(13) = {0.15, 2.5, 0.1, h};
Point(14) = {0.15, 2.5, 0.9, h};
Point(15) = {-0.15, 2.5, 0.1, h};
Point(16) = {-0.15,2.5, 0.9, h};

Circle(1) = {8, 2, 6};
Circle(2) = {12, 4, 10};
Circle(3) = {11, 3, 9};
Circle(4) = {7, 1, 5};
Line(5) = {6, 10};
Line(6) = {8, 12};
Line(7) = {10, 9};
Line(8) = {9, 5};
Line(9) = {7, 11};
Line(10) = {11, 12};
Line(11) = {10, 14};
Line(12) = {12, 16};
Line(13) = {11, 15};
Line(14) = {9, 13};
Line(15) = {13, 14};
Line(16) = {16, 15};
Line(17) = {15, 13};
Line(18) = {16, 14};


//BOTTOM BLADE

//Bottom part
Point(17) = {0.15, -0.476969, 0, h};
Point(18) = {0.15, -0.476969, 1, h};
Point(19) = {-0.15, -0.476969, 0, h};
Point(20) = {-0.15, -0.476969, 1, h};

Point(21) = {0.15, -0.476969, 0.1, h};
Point(22) = {0.15, -0.476969, 0.9, h};
Point(23) = {-0.15, -0.476969, 0.1, h};
Point(24) = {-0.15, -0.476969, 0.9, h};

//Top part

Point(25) = {0.15, -2.5, 0.1, h};
Point(26) = {0.15, -2.5, 0.9, h};
Point(27) = {-0.15, -2.5, 0.1, h};
Point(28) = {-0.15,-2.5, 0.9, h};

Circle(19) = {20, 2, 18};
Circle(20) = {24, 4, 22};
Circle(21) = {23, 3, 21};
Circle(22) = {19, 1, 17};
Line(23) = {18, 22};
Line(24) = {22, 21};
Line(25) = {21, 17};
Line(26) = {20, 24};
Line(27) = {24, 23};
Line(28) = {23, 19};
Line(29) = {22, 26};
Line(30) = {24, 28};
Line(31) = {21, 25};
Line(32) = {23, 27};
Line(33) = {27, 28};
Line(34) = {28, 26};
Line(35) = {26, 25};
Line(36) = {25, 27};

//LEFT BLADE

//Bottom part
Point(29) = {-0.476969, 0.15, 0, h};
Point(30) = {-0.476969, 0.15, 1, h};
Point(31) = {-0.476969,-0.15 , 0, h};
Point(32) = {-0.476969,-0.15 , 1, h};

Point(33) = {-0.476969, 0.15, 0.1, h};
Point(34) = {-0.476969, 0.15, 0.9, h};
Point(35) = {-0.476969,-0.15 , 0.1, h};
Point(36) = {-0.476969,-0.15 , 0.9, h};

//Top part

Point(37) = {-2.5, 0.15, 0.1, h};
Point(38) = {-2.5, 0.15, 0.9, h};
Point(39) = {-2.5,-0.15 , 0.1, h};
Point(40) = {-2.5,-0.15 , 0.9, h};

Circle(37) = {32, 2, 30};
Circle(38) = {36, 4, 34};
Circle(39) = {35, 3, 33};
Circle(40) = {31, 1, 29};
Line(41) = {29, 33};
Line(42) = {33, 34};
Line(43) = {34, 30};
Line(44) = {32, 36};
Line(45) = {36, 35};
Line(46) = {35, 31};
Line(47) = {36, 40};
Line(48) = {38, 34};
Line(49) = {33, 37};
Line(50) = {35, 39};
Line(51) = {37, 38};
Line(52) = {38, 40};
Line(53) = {40, 39};
Line(54) = {39, 37};

//RIGHT BLADE

//Bottom part
Point(41) = {0.476969, 0.15, 0, h};
Point(42) = {0.476969, 0.15, 1, h};
Point(43) = {0.476969,-0.15 , 0, h};
Point(44) = {0.476969,-0.15 , 1, h};

Point(45) = {0.476969, 0.15, 0.1, h};
Point(46) = {0.476969, 0.15, 0.9, h};
Point(47) = {0.476969,-0.15 , 0.1, h};
Point(48) = {0.476969,-0.15 , 0.9, h};

//Top part

Point(49) = {2.5, 0.15, 0.1, h};
Point(50) = {2.5, 0.15, 0.9, h};
Point(51) = {2.5,-0.15 , 0.1, h};
Point(52) = {2.5,-0.15 , 0.9, h};
//+
Circle(55) = {44, 2, 42};
//+
Circle(56) = {48, 4, 46};
//+
Circle(57) = {45, 3, 47};
//+
Circle(58) = {41, 1, 43};
//+
Line(59) = {41, 45};
//+
Line(60) = {45, 46};
//+
Line(61) = {46, 42};
//+
Line(62) = {43, 47};
//+
Line(63) = {47, 48};
//+
Line(64) = {48, 44};
//+
Line(65) = {50, 46};
//+
Line(66) = {52, 48};
//+
Line(67) = {51, 52};
//+
Line(68) = {49, 50};
//+
Line(69) = {50, 52};
//+
Line(70) = {51, 49};
//+
Line(71) = {49, 45};
//+
Line(72) = {51, 47};
//+
Circle(73) = {32, 2, 30};
//+
Circle(74) = {20, 2, 32};
//+
Circle(75) = {2, 30, 8};
//+
Circle(75) = {30, 2, 8};
//+
Circle(76) = {6, 2, 42};
//+
Circle(77) = {44, 2, 18};
//+
Circle(78) = {41, 1, 5};
//+
Circle(79) = {7, 1, 29};
//+
Circle(80) = {31, 1, 19};
//+
Circle(81) = {17, 1, 43};
//+
Delete {
  Curve{8}; Curve{9}; Curve{23}; Curve{5}; Curve{6}; Curve{59}; Curve{61}; Curve{62}; Curve{64}; Curve{25}; Curve{28}; Curve{46}; Curve{44}; Curve{26}; Curve{43}; Curve{41}; 
}
//+
Line(82) = {11, 3};
//+
Delete {
  Curve{82}; 
}
//+
Line(82) = {11, 7};
//+
Line(83) = {9, 5};
//+
Line(84) = {45, 41};
//+
Line(85) = {47, 43};
//+
Line(86) = {21, 17};
//+
Line(87) = {23, 19};
//+
Line(88) = {31, 35};
//+
Line(89) = {29, 33};
//+
Line(90) = {34, 30};
//+
Line(91) = {36, 32};
//+
Line(92) = {20, 24};
//+
Line(93) = {22, 18};
//+
Line(94) = {44, 48};
//+
Line(95) = {42, 46};
//+
Line(96) = {6, 10};
//+
Line(97) = {8, 12};
//+
Curve Loop(1) = {33, 34, 35, 36};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {31, -35, -29, 24};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {30, 34, -29, -20};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {32, 33, -30, 27};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {32, -36, -31, -21};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {50, 54, -49, -39};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {47, 53, -50, -45};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {52, 53, 54, 51};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {48, -42, 49, 51};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {47, -52, 48, -38};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {13, 17, -14, -3};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {16, 17, 15, -18};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {12, 16, -13, 10};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {14, 15, -11, 7};
//+
Plane Surface(14) = {14};
//+
Curve Loop(15) = {11, -18, -12, 2};
//+
Plane Surface(15) = {15};
//+
Curve Loop(16) = {65, -56, -66, -69};
//+
Plane Surface(16) = {16};
//+
Curve Loop(17) = {68, 69, -67, 70};
//+
Plane Surface(17) = {17};
//+
Curve Loop(18) = {67, 66, -63, -72};
//+
Plane Surface(18) = {18};
//+
Curve Loop(19) = {57, -72, 70, 71};
//+
Plane Surface(19) = {19};
//+
Curve Loop(20) = {65, -60, -71, 68};
//+
Plane Surface(20) = {20};
//+
Curve Loop(21) = {78, -4, 79, -40, 80, 22, 81, -58};
//+
Plane Surface(21) = {21};
//+
Curve Loop(22) = {77, -19, 74, 37, 75, 1, 76, -55};
//+
Plane Surface(22) = {22};
//+
Curve Loop(23) = {57, 63, 56, -60};
//+
Surface(23) = {23};
//+
Curve Loop(25) = {95, -56, -94, 55};
//+
Surface(24) = {25};
//+
Curve Loop(27) = {58, -85, -57, 84};
//+
Surface(25) = {27};
//+
Curve Loop(29) = {78, -83, -7, -96, 76, 95, -60, 84};
//+
Surface(26) = {29};
//+
Curve Loop(31) = {96, -2, -97, 1};
//+
Surface(27) = {31};
//+
Curve Loop(33) = {83, -4, -82, 3};
//+
Surface(28) = {33};
//+
Curve Loop(35) = {79, 89, 42, 90, 75, 97, -10, 82};
//+
Surface(29) = {35};
//+
Curve Loop(37) = {90, -37, -91, 38};
//+
Surface(30) = {37};
//+
Curve Loop(39) = {89, -39, -88, 40};
//+
Surface(31) = {39};
//+
Curve Loop(41) = {80, -87, -27, -92, 74, -91, 45, -88};
//+
Surface(32) = {41};
//+
Curve Loop(43) = {92, 20, 93, -19};
//+
Surface(33) = {43};
//+
Curve Loop(45) = {22, -86, -21, 87};
//+
Surface(34) = {45};
//+
Curve Loop(47) = {42, -38, 45, 39};
//+
Surface(35) = {47};
//+
Curve Loop(49) = {27, 21, -24, -20};
//+
Surface(36) = {49};
//+
Curve Loop(51) = {86, 81, -85, 63, -94, 77, -93, 24};
//+
Surface(37) = {51};
//+
Curve Loop(53) = {3, -7, -2, -10};
//+
Surface(38) = {53};
//+
Surface Loop(1) = {3, 4, 5, 1, 2, 36};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {7, 10, 8, 6, 9, 35};
//+
Volume(2) = {2};
//+
Surface Loop(3) = {15, 14, 11, 13, 12, 38};
//+
Volume(3) = {3};
//+
Surface Loop(4) = {19, 18, 17, 20, 16, 23};
//+
Volume(4) = {4};
//+
Surface Loop(5) = {26, 21, 28, 29, 31, 32, 34, 37, 33, 22, 30, 27, 24, 25, 23, 38, 35, 36};
//+
Volume(5) = {5};
//+
Physical Volume("Left blade", 98) = {2};
//+
Physical Volume("Top blade", 99) = {3};
//+
Physical Volume("Right blade", 100) = {4};
//+
Physical Volume("Bottom blade", 101) = {1};
//+
Physical Volume("Center", 102) = {5};
