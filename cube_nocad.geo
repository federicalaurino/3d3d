// Box in a box
SetFactory("Built-in");

DefineConstant[
wx = {0.25, Name "inner x half width"}
wy = {0.25, Name "inner y half widht"}
Wx = {0.5, Name "outer x half width"}
Wy = {0.5, Name "outer y half width"}
H = {1, Name "height"}
];


size = 0.125;
Point(1) = {-wx, -wy, 0, size};
Point(2) = {wx, -wy, 0, size};
Point(3) = {wx, wy, 0, size};
Point(4) = {-wx, wy, 0, size};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};  // inner

Point(5) = {-Wx, -Wy, 0, size};
Point(6) = {Wx, -Wy, 0, size};
Point(7) = {Wx, Wy, 0, size};
Point(8) = {-Wx, Wy, 0, size};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {1, 2};

surfs = {1, 2};
Extrude {0, 0, H} {
  Surface{surfs[]}; 
//for structured meshes
//Layers{30}; 
}
Physical Volume(70) = {1};
Physical Volume(7) = {2};

// Outer bdries
Physical Surface(1) = {67};
Physical Surface(2) = {63};
Physical Surface(3) = {59};
Physical Surface(4) = {71};
Physical Surface(5) = {72};
Physical Surface(6) = {2};

// // Inner bdries
Physical Surface(50) = {30};
Physical Surface(60) = {1};
Physical Surface(10) = {21};
Physical Surface(20) = {25};
Physical Surface(30) = {29};
Physical Surface(40) = {17};
//
