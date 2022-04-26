// Gmsh project created on Sat Mar 19 17:07:45 2022

Mesh.MshFileVersion = 2.2;
Mesh.SubdivisionAlgorithm = 2;

h=0.7;


//CENTERS
Point(1) = {0,0,0,h};
Point(2) = {1,0,0,h};


//TOP BLADE

//Bottom part
Point(3) = {0.05, 0.476969, 0.15, h};
Point(4) = {0.95, 0.476969, 0.15, h};
Point(5) = {0.05, 0.476969, -0.15, h};
Point(6) = {0.95, 0.476969, -0.15, h};

//Top part
Point(7) = {0.95,2.5,0.15,h};
Point(8) = {0.95,2.5,-0.15,h};
Point(9) = {0.05,2.5,0.15,h};
Point(10) = {0.05,2.5,-0.15,h};

Line(1) = {3,9};
Line(2) = {4,7};
Line(3) = {5,10};
Line(4) = {6,8};
Line(5) = {3,4};
Line(6) = {5,6};
Line(7) = {7,9};
Line(8) = {8,10};
Line(9) = {9,3};
Line(10) = {7,8};
Line(11) = {9,10};
Circle(12) = {5,1,3};
Circle(13) = {6,2,4};

Curve Loop(1) = {3, -11, -1, -12};
Surface(1) = {1};
Curve Loop(2) = {5, 2, 7, -1};
Surface(2) = {2};
Curve Loop(3) = {3, -8, -4, -6};
Surface(3) = {3};
Curve Loop(4) = {11, -8, -10, 7};
Surface(4) = {4};
Curve Loop(5) = {2, 10, -4, 13};
Surface(5) = {5};
Curve Loop(6) = {12, 5, -13, -6};
Surface(6) = {6};


Surface Loop(1) = {4, 1, 3, 5, 2, 6};
Volume(1) = {1};


//Bottom blade

//Top part
Point(11) = {0, -0.476969, 0.15, h};
Point(12) = {1, -0.476969, 0.15, h};
Point(13) = {0, -0.476969, -0.15, h};
Point(14) = {1, -0.476969, -0.15, h};

//Bottom part
Point(15) = {1,-2.5,0.15,h};
Point(16) = {1,-2.5,-0.15,h};
Point(17) = {0,-2.5,0.15,h};
Point(18) = {0,-2.5,-0.15,h};

Line(27) = {15, 16};
Line(28) = {17, 18};
Line(29) = {18, 16};
Line(30) = {15, 17};
Line(31) = {15, 12};
Line(32) = {16, 14};
Line(33) = {18, 13};
Line(34) = {17, 11};
Line(35) = {13, 14};
Line(36) = {11, 12};
Circle(37) = {12, 2, 14};
Circle(38) = {13, 1, 11};

Curve Loop(13) = {33, 35, -32, -29};
Plane Surface(13) = {13};
Curve Loop(14) = {30, 28, 29, -27};
Plane Surface(14) = {14};
Curve Loop(15) = {34, 36, -31, 30};
Plane Surface(15) = {15};
Curve Loop(16) = {31, 37, -32, -27};
Plane Surface(16) = {16};
Curve Loop(17) = {28, 33, 38, -34};
Plane Surface(17) = {17};
Curve Loop(18) = {36, 37, -35, 38};
Surface(18) = {18};

Surface Loop(3) = {13, 17, 14, 15, 18, 16};
Volume(3) = {3};



//Left blade



Point(19) = {0, 0.15, 0.476969, h};
Point(20) = {1, 0.15, 0.476969, h};
Point(21) = {0, -0.15, 0.476969, h};
Point(22) = {1, -0.15, 0.476969, h};

Point(23) = {1,0.15,2.5,h};
Point(24) = {1,-0.15,2.5,h};
Point(25) = {0,0.15,2.5,h};
Point(26) = {0,-0.15,2.5,h};

Line(14) = {20, 23};
Line(15) = {22, 24};
Line(16) = {19, 25};
Line(17) = {21, 26};
Line(18) = {19, 19};
Line(19) = {20, 19};
Line(20) = {21, 22};
Line(21) = {24, 23};
Line(22) = {26, 25};
Line(23) = {25, 23};
Line(24) = {26, 24};
Circle(25) = {22, 2, 20};
Circle(26) = {21, 1, 19};

Curve Loop(7) = {16, 23, -14, 19};
Plane Surface(7) = {7};
Curve Loop(8) = {23, -21, -24, 22};
Plane Surface(8) = {8};
Curve Loop(9) = {24, -15, -20, 17};
Plane Surface(9) = {9};
Curve Loop(10) = {14, -21, -15, 25};
Surface(10) = {10};
Curve Loop(11) = {26, 16, -22, -17};
Surface(11) = {11};
Curve Loop(12) = {20, 25, 19, -26};
Surface(12) = {12};

Surface Loop(2) = {7, 11, 12, 9, 8, 10};
Volume(2) = {2};


//Right blade



Point(27) = {0, 0.15, -0.476969, h};
Point(28) = {1, 0.15, -0.476969, h};
Point(29) = {0, -0.15, -0.476969, h};
Point(30) = {1, -0.15, -0.476969, h};

Point(31) = {1,0.15,-2.5,h};
Point(32) = {1,-0.15,-2.5,h};
Point(33) = {0,0.15,-2.5,h};
Point(34) = {0,-0.15,-2.5,h};

Line(39) = {28, 31};
Line(40) = {30, 32};
Line(41) = {32, 31};
Line(42) = {31, 33};
Line(43) = {33, 34};
Line(44) = {34, 32};
Line(45) = {34, 29};
Line(46) = {33, 27};
Line(47) = {27, 28};
Line(48) = {29, 30};
Circle(49) = {28, 2, 30};
Circle(50) = {27, 1, 29};

Curve Loop(19) = {46, 47, 39, 42};
Plane Surface(19) = {19};
Curve Loop(20) = {42, 43, 44, 41};
Plane Surface(20) = {20};
Curve Loop(21) = {39, -41, -40, -49};
Plane Surface(21) = {21};
Curve Loop(22) = {46, 50, -45, -43};
Plane Surface(22) = {22};
Curve Loop(23) = {47, 49, -48, -50};
Surface(23) = {23};
Curve Loop(24) = {45, 48, 40, -44};
Plane Surface(24) = {24};

Surface Loop(4) = {19, 22, 23, 21, 20, 24};
Volume(4) = {4};

//Heart

Circle(51) = {4, 2, 20};
Circle(52) = {6, 2, 28};
Circle(53) = {30, 2, 14};
Circle(54) = {12, 2, 22};
Circle(55) = {5, 1, 27};
Circle(56) = {13, 1, 29};
Circle(57) = {11, 1, 21};
Circle(58) = {19, 1, 3};

Curve Loop(25) = {52, 49, 53, -37, 54, 25, -51, -13};
Plane Surface(25) = {25};
Curve Loop(26) = {55, 50, -56, 38, 57, 26, 58, -12};
Plane Surface(26) = {26};
Curve Loop(27) = {6, 52, -47, -55};
Surface(27) = {27};
Curve Loop(28) = {58, 5, 51, 19};
Surface(28) = {28};
Curve Loop(29) = {48, 53, -35, 56};
Surface(29) = {29};
Curve Loop(30) = {57, 20, -54, -36};
Surface(30) = {30};

Surface Loop(5) = {25, 27, 26, 29, 30, 28, 23, 12, 18, 6};
Volume(5) = {5};


//Physical volumes

Physical Volume("Left blade", 59) = {2};
Physical Volume("Top blade", 60) = {1};
Physical Volume("Right blade", 61) = {4};
Physical Volume("Bottom blade", 62) = {3};
Physical Volume("Heart", 63) = {5};
