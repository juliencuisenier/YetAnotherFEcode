h=0.5;


Point(1) = {0,0,0,h};
Point(2) = {0.33,0,0,h};
Point(3) = {0.66,0,0,h};
Point(4) = {1,0,0,h};
Point(5) = {0,0.5,0,h};
Point(6) = {0,0.5,0.3,h};
Point(7) = {0.33,0.5,0.3,h};
Point(8) = {0.33,0.5,0,h};
Point(9) = {0.66,0.5,0.3,h};
Point(10) = {0.66,0.5,0,h};
Point(11) = {0.66,0,0.3,h};
Point(12) = {1,0.5,0.3,h};
Point(13) = {1,0,0.3,h};
Point(14) = {1,0.5,0,h};
Point(15) = {0,0,0.3,h};
Point(16) = {0.33,0,0.3,h};//+
Line(1) = {5, 8};
//+
Line(2) = {8, 7};
//+
Line(3) = {7, 6};
//+
Line(4) = {6, 5};
//+
Line(5) = {6, 15};
//+
Line(6) = {7, 16};
//+
Line(7) = {16, 15};
//+
Line(8) = {15, 1};
//+
Line(9) = {1, 2};
//+
Line(10) = {2, 8};
//+
Line(11) = {16, 2};
//+
Line(12) = {1, 5};
//+
Line(13) = {7, 9};
//+
Line(14) = {8, 10};
//+
Line(15) = {10, 3};
//+
Line(16) = {3, 11};
//+
Line(17) = {11, 16};
//+
Line(18) = {11, 9};
//+
Line(19) = {9, 10};
//+
Line(20) = {10, 14};
//+
Line(21) = {14, 12};
//+
Line(22) = {12, 9};
//+
Line(23) = {4, 14};
//+
Line(24) = {12, 13};
//+
Line(25) = {11, 13};
//+
Line(26) = {13, 4};
//+
Line(27) = {4, 3};
//+
Line(28) = {3, 2};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {12, -4, 5, 8};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {7, -5, -3, 6};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {14, 15, 28, 10};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {14, -19, -13, -2};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {18, 19, 15, 16};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {17, 11, -28, 16};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {27, 16, 25, 26};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {20, -23, 27, -15};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {21, 22, 19, 20};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {24, 26, 23, 21};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {6, -17, 18, -13};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {22, -18, 25, -24};
//+
Plane Surface(13) = {13};
//+
Surface Loop(1) = {9, 10, 11, 13, 8, 6};
//+
Volume(1) = {1};
//+
Curve Loop(14) = {10, 2, 6, 11};
//+
Plane Surface(14) = {14};
//+
Surface Loop(2) = {5, 4, 7, 12, 14, 6};
//+
Volume(2) = {2};
//+
Curve Loop(15) = {1, -10, -9, 12};
//+
Plane Surface(15) = {15};
//+
Curve Loop(16) = {8, 9, -11, 7};
//+
Plane Surface(16) = {16};
//+
Surface Loop(3) = {15, 1, 2, 3, 16, 14};
//+
Volume(3) = {3};
