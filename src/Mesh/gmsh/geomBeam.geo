Mesh.MshFileVersion = 2.2;

// definition du pas du maillage
h = 0.1;

// definition des points (en 3D, raison pour laquelle il y a un 0 en z)
Point(1) = {0, 0, 0, h};
Point(2) = {0.2, 0, 0, h};
Point(3) = {0, 1, 0, h};
Point(4) = {0.2, 1, 0, h};
Point(5) = {0, 0, 0.3, h};
Point(6) = {0.2, 0, 0.3, h};
Point(7) = {0, 1, 0.3, h};
Point(8) = {0.2, 1, 0.3, h};
Point(9) = {0, 0.5, 0, h};
Point(10) = {0, 0.5, 0.3, h};
Point(11) = {0.2, 0.5, 0, h};
Point(12) = {0.2, 0.5, 0.3, h};

// definition des segments qui relient les points
Line(1) = {1,2};
Line(2) = {1,5};
Line(3) = {2,6};
Line(4) = {5,6};
Line(5) = {5,10};
Line(6) = {6,12};
Line(7) = {1,9};
Line(8) = {2,11};
Line(9) = {12,11};
Line(10) = {12,10};
Line(11) = {10,9};
Line(12) = {11,9};
Line(13) = {4,11};
Line(14) = {12,8};
Line(15) = {9,3};
Line(16) = {10,7};
Line(17) = {7,3};
Line(18) = {7,8};
Line(19) = {3,4};
Line(20) = {4,8};

//+
Curve Loop(1) = {2, 4, -3, -1};
//+
Surface(1) = {1};
//+
Curve Loop(2) = {5, -10, -6, -4};
//+
Surface(2) = {2};
//+
Curve Loop(3) = {11, -12, -9, 10};
//+
Surface(3) = {3};
//+
Curve Loop(4) = {7, -12, -8, -1};
//+
Surface(4) = {4};
//+
Curve Loop(5) = {15, 19, 13, 12};
//+
Surface(5) = {5};
//+
Curve Loop(6) = {20, -18, 17, 19};
//+
Surface(6) = {6};
//+
Curve Loop(7) = {16, 18, -14, 10};
//+
Surface(7) = {7};
//+
Curve Loop(8) = {16, 17, -15, -11};
//+
Surface(8) = {8};
//+
Curve Loop(9) = {9, -13, 20, -14};
//+
Surface(9) = {9};
//+
Curve Loop(10) = {8, -9, -6, -3};
//+
Surface(10) = {10};
//+
Curve Loop(11) = {5, 11, -7, 2};
//+
Surface(11) = {11};
//+
Surface Loop(1) = {2, 11, 4, 10, 1, 7, 8, 6, 9, 5};
//+
Surface Loop(2) = {2, 11, 4, 10, 1, 3};
//+
Volume(1) = {2};
//+
Surface Loop(3) = {7, 8, 6, 9, 5, 3};
//+
Volume(2) = {3};
//+
Physical Point("Beam1", 21) = {7, 8, 12, 4, 3, 9, 11, 10};
//+
Physical Point("Beam2", 22) = {12, 6, 5, 10, 9, 11, 2, 1};
//+
Physical Surface("Beam1", 23) = {7, 6, 9, 8, 5, 3};
//+
Physical Surface("Beam2", 24) = {2, 10, 11, 3, 1, 4};
//+
Physical Volume("Beam1", 25) = {2};
//+
Physical Volume("Beam2", 26) = {1};
