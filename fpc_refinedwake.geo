// Gmsh project
SetFactory("OpenCASCADE");

// ---------------- Geometry ----------------
Point(1) = {0,    0,    0, 1.0};
Point(2) = {0,    0.41, 0, 1.0};
Point(3) = {2.2,  0.41, 0, 1.0};
Point(4) = {2.2,  0,    0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
// (fixed) remove zero-length edge; keep the outer boundary close
Line(4) = {4, 1};

// Cylinder (hole) center (0.2, 0.2), radius 0.05
Circle(5) = {0.2, 0.2, 0, 0.05, 0, 2*Pi};

// Physical groups (keep your IDs)
Physical Curve(4) = {1}; // left
Physical Curve(1) = {2}; // top
Physical Curve(2) = {3}; // right
Physical Curve(3) = {4}; // bottom
Physical Curve(5) = {5}; // cylinder wall

Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5};
Plane Surface(1) = {1, 2};

Physical Surface(6) = {1};
// ---------------- targeted sizing: cylinder + short wake ----------------
xc = 0.2; yc = 0.2; R = 0.05;        // cylinder center and radius
Lx = 2.2; H = 0.41;                  // domain dimensions
D  = 2*R;

// --- base mesh size (uniform, coarse elsewhere)
lcBase = 0.125;                       // coarse size away from cylinder/wake
Field[1] = MathEval;                 // constant background size
Field[1].F = Sprintf("%g", lcBase);

// --- (A) cylinder ring refinement (Distance + Threshold)
lcCyl   = 0.05;                     // finest size at the wall
bandOut = 0.20;                      // how far the refinement ring extends

Field[10] = Distance;
Field[10].CurvesList = {5};          // the cylinder circular curve id
Field[10].Sampling   = 100;

Field[11] = Threshold;
Field[11].InField = 10;
Field[11].SizeMin = lcCyl;           // inside band
Field[11].SizeMax = lcBase;          // smoothly back to base size
Field[11].DistMin = R;               // from cylinder boundary
Field[11].DistMax = R + bandOut;     // to this distance

// --- (B) wake refinement: width ~ (widthFactor * D), length up to mid-channel
widthFactor   = 1.3;                 // >1 → a bit wider than the diameter
wakeHalf      = 1*widthFactor*D;   // half thickness of wake strip
lcWake        = 0.05;               // target size inside the strip

Field[20] = Box;
Field[20].VIn  = lcWake;
Field[20].VOut = lcBase;
Field[20].XMin = xc + 1.05*R;        // start just after cylinder
Field[20].XMax = Lx/3;               // stop at channel midpoint
Field[20].YMin = yc - wakeHalf;      // narrow band around the centerline
Field[20].YMax = yc + wakeHalf;

// --- Combine fields: smallest size wins
Field[100] = Min;
Field[100].FieldsList = {1, 11, 20};
Background Field = 100;

// Global bounds (optional, keeps mesh generator in check)
Mesh.CharacteristicLengthMin = lcCyl;
Mesh.CharacteristicLengthMax = lcBase;

