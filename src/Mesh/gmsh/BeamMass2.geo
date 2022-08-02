// Gmsh project created on Wed Jun 29 17:23:49 2022
SetFactory("OpenCASCADE");



Mesh.MshFileVersion = 2.2;
Mesh.SubdivisionAlgorithm = 2;

h = 1;


//Beam

Point(1) = {0,0.15,0.25,h};
Point(2) = {0,-0.15,-0.25,h};
Point(3) = {0,0.15,-0.25,h};
Point(4) = {0,-0.15,0.25,h};
Point(5) = {5,0.15,0.25,h};
Point(6) = {5,-0.15,-0.25,h};
Point(7) = {5,0.15,-0.25,h};
Point(8) = {5,-0.15,0.25,h};

Line(1) = {3, 7};
Line(2) = {7, 6};
Line(3) = {6, 8};
Line(4) = {8, 5};
Line(5) = {5, 1};
Line(6) = {4, 8};
Line(7) = {5, 7};
Line(8) = {3, 1};
Line(9) = {1, 4};
Line(10) = {4, 2};
Line(11) = {2, 3};
Line(12) = {2, 6};

//Mass

Point(9) = {5, 0.5,0.5,h};
Point(10) = {5, -0.5,0.5,h};
Point(11) = {5, 0.5,-0.5,h};
Point(12) = {5, -0.5,-0.5,h};
Point(13) = {6, 0.5,0.5,h};
Point(14) = {6, -0.5,0.5,h};
Point(15) = {6, 0.5,-0.5,h};
Point(16) = {6, -0.5,-0.5,h};//+
Line(13) = {9, 11};
//+
Line(14) = {11, 15};
//+
Line(15) = {15, 16};
//+
Line(16) = {16, 14};
//+
Line(17) = {14, 10};
//+
Line(18) = {10, 12};
//+
Line(19) = {12, 11};
//+
Line(20) = {9, 13};
//+
Line(21) = {13, 14};
//+
Line(22) = {10, 9};
//+
Line(23) = {13, 15};
//+
Curve Loop(1) = {1, -7, 5, -8};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, 9, 6, 4};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {2, 3, 4, 7};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {1, 2, -12, 11};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {8, 9, 10, 11};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {13, 14, -23, -20};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {15, 16, -21, 23};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {17, 22, 20, 21};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {13, -19, -18, 22};
//+
Plane Surface(9) = {9};
//+
Line(24) = {12, 16};
//+
Curve Loop(10) = {14, 15, -24, 19};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {6, -3, -12, -10};
//+
Plane Surface(11) = {11};


//+
Curve Loop(12) = {18, 24, 16, 17};
//+
Plane Surface(12) = {12};
//+
Recursive Delete {
  Surface{9}; 
}
//+
Recursive Delete {
  Surface{3}; 
}
//+
Curve Loop(13) = {13, -19, -18, 22};
//+
Curve Loop(14) = {4, 7, 2, 3};
//+
Plane Surface(13) = {13, 14};
//+
Curve Loop(15) = {7, 2, 3, 4};
//+
Surface Loop(1) = {1, 4, 13, 6, 10, 7, 12, 8, 11, 2, 5};
//+
Volume(1) = {1};
//+
Physical Volume("system", 25) = {1};
