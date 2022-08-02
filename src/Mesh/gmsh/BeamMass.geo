// Gmsh project created on Wed Jun 29 13:44:08 2022
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {-1, -0.1, 0, 5, 0.1, 0};
//+
Extrude {0, 0, 0.3} {
  Curve{3}; Curve{4}; Curve{1}; Curve{2}; 
}

Mesh.MshFileVersion = 2.2;
Mesh.SubdivisionAlgorithm = 2;

h=2;

Point(9) = {4, 0.45, 0.65,h};
Point(10) = {4,0.45,-0.35,h};
Point(11) = {4, -0.55,0.65,h};
Point(12) = {4,-0.55,-0.35,h};
//+
Line(13) = {9, 10};
//+
Line(14) = {10, 12};
//+
Line(15) = {12, 11};
//+
Line(16) = {11, 9};
//+
Extrude {1, 0, 0} {
  Curve{13}; Curve{14}; Curve{15}; Curve{16}; 
}


//+
Curve Loop(10) = {19, 21, 23, 24};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {13, 14, 15, 16};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {7, 9, 11, 12};
//+
Plane Surface(12) = {12};
//+
Surface Loop(1) = {4, 3, 2, 5, 12, 1};
//+
Surface Loop(2) = {6, 9, 8, 7, 10, 11};
//+
Volume(1) = {1, 2};
