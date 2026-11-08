(* vector_field_spherical_coordinate_transform.wl
   Matrix pieces used by python/ads/create_vector_field_rotation_vars.py to
   rotate a radiated vector field E at a fixed observation point r:

     E'(r) = R . E(R^-1 . r)

   i.e. sample the unrotated pattern at the body-frame direction R^-1.r and
   rotate the field vector by R. Writes to tables/:
     Cart2Sph_ADS_Equations.txt  ADS Eqn text: Cart2Sph_ij(theta,phi) at the
                                 observation direction, and the generic
                                 Cartesian-ratio template M_ij(x,y,z,r,rho) with
                                 r_of_xyz, rho_of_xyz, theta_of_xyz, phi_of_xyz.
   In ADS the primed (body-frame) quantities are named with a _prime suffix. *)

ClearAll[Cart2SphSymbolic, Sph2CartSymbolic, Cart2SphCartesianSymbolic,
  rDef, rhoDef, GetCart2Sphij, GetSph2Cartij, GetCart2SphCartesianij,
  AdsString, BuildMatrixIndexedEntries, BuildADSEquationFile,
  BuildPreviewGrid];

(* Run from the mathematica/ folder; generated tables are written to tables/. *)
SetDirectory[NotebookDirectory[]];
If[!DirectoryQ["tables"], CreateDirectory["tables"]];
Print["Working directory set to: ", Directory[]];

(* ============================================================ *)
(* Core symbolic definition -- angle form                        *)
(* ============================================================ *)

(* Cart2Sph is the orthogonal change-of-basis matrix whose ROWS are
   the spherical unit vectors rhat, thetahat, phihat expressed in
   Cartesian components (physics convention: theta = polar angle from
   +z, phi = azimuthal angle in the xy-plane):

     rhat     = ( sin(theta)cos(phi),  sin(theta)sin(phi),  cos(theta))
     thetahat = ( cos(theta)cos(phi),  cos(theta)sin(phi), -sin(theta))
     phihat   = (-sin(phi),            cos(phi),            0        )

   so Cart2Sph.{Vx,Vy,Vz} gives {Vr,Vtheta,Vphi}. Sph2Cart, the reverse
   map, is its transpose -- valid since Cart2Sph is orthogonal. *)
Cart2SphSymbolic = {
   {Sin[theta]*Cos[phi], Sin[theta]*Sin[phi], Cos[theta]},
   {Cos[theta]*Cos[phi], Cos[theta]*Sin[phi], -Sin[theta]},
   {-Sin[phi], Cos[phi], 0}
   };

Sph2CartSymbolic = Transpose[Cart2SphSymbolic];

GetCart2Sphij[i_, j_] := Cart2SphSymbolic[[i, j]];
GetSph2Cartij[i_, j_] := Sph2CartSymbolic[[i, j]];

(* ============================================================ *)
(* Core symbolic definition -- CARTESIAN form of Cart2Sph          *)
(* ============================================================ *)

(* Sph2Cart is needed at the body-frame direction (theta',phi') of
   (x',y',z') = R^T.(x,y,z), where (x,y,z) = r*rhat(theta,phi) is the fixed
   observation position. To avoid inverse trig, the trig functions are
   written as ratios of x,y,z (r = sqrt(x^2+y^2+z^2), rho = sqrt(x^2+y^2)):

     sin(theta) = rho/r     cos(theta) = z/r
     sin(phi)   = y/rho     cos(phi)   = x/rho

   Substituting these into Cart2SphSymbolic and simplifying gives:

     rhat     = ( x/r,          y/r,          z/r   )
     thetahat = ( x*z/(r*rho),  y*z/(r*rho), -rho/r  )
     phihat   = (-y/rho,        x/rho,        0      )

   This generic template is exported as M_ij; the Python script uses it
   transposed (M_ji) with x,y,z,r,rho -> x_prime,...,rho_prime to build
   Sph2Cart at the body-frame point (Sph2Cart(p) = Cart2Sph(p)^T).
   r and rho are kept as free symbols (rDef, rhoDef) so each square root is
   computed once in ADS. *)
rDef = Sqrt[x^2 + y^2 + z^2];
rhoDef = Sqrt[x^2 + y^2];

Cart2SphCartesianSymbolic = {
   {x/r, y/r, z/r},
   {x*z/(r*rho), y*z/(r*rho), -rho/r},
   {-y/rho, x/rho, 0}
   };

GetCart2SphCartesianij[i_, j_] := Cart2SphCartesianSymbolic[[i, j]];

(* ============================================================ *)
(* Sanity checks                                                 *)
(* ============================================================ *)

SeedRandom[7];
thetaVal = RandomReal[{0, Pi}];
phiVal = RandomReal[{0, 2 Pi}];
subsVal = {theta -> thetaVal, phi -> phiVal};

(* 1) Orthogonality: Cart2Sph . Cart2Sph^T should be the identity,
   since its rows are an orthonormal triad. *)
orthoCheck = N[Cart2SphSymbolic /. subsVal] .
   Transpose[N[Cart2SphSymbolic /. subsVal]];
orthoMaxDiff = Max[Abs[Flatten[orthoCheck - IdentityMatrix[3]]]];
Print["Orthogonality check -- max |Cart2Sph.Cart2Sph^T - I| = ",
  orthoMaxDiff, " (should be ~0)"];

(* 2) Handedness: Det[Cart2Sph] should be +1. *)
detVal = Det[N[Cart2SphSymbolic /. subsVal]];
Print["Handedness check -- Det[Cart2Sph] = ", detVal, " (should be +1)"];

(* 3) Reference point: at theta=Pi/2, phi=0, rhat=(1,0,0),
   thetahat=(0,0,-1), phihat=(0,1,0). *)
refCheck = N[Cart2SphSymbolic /. {theta -> Pi/2, phi -> 0}];
Print["Reference-point check at theta=Pi/2,phi=0 (rows = rhat,thetahat,phihat): ",
  refCheck, " (expected {{1,0,0},{0,0,-1},{0,1,0}})"];

(* 4) Cartesian-form check: Cart2SphCartesianSymbolic, evaluated at the
   (x,y,z) = rhat(theta,phi) corresponding to a random test point,
   should agree with Cart2SphSymbolic(theta,phi) directly -- confirming
   the trig-identity rewrite above is correct, independent of the
   angle-based formula. *)
xVal = Sin[thetaVal]*Cos[phiVal];
yVal = Sin[thetaVal]*Sin[phiVal];
zVal = Cos[thetaVal];
rVal = rDef /. {x -> xVal, y -> yVal, z -> zVal};
rhoVal = rhoDef /. {x -> xVal, y -> yVal, z -> zVal};
cartCheck = N[Cart2SphCartesianSymbolic /. {x -> xVal, y -> yVal,
    z -> zVal, r -> rVal, rho -> rhoVal}];
cartMaxDiff = Max[Abs[Flatten[cartCheck - N[Cart2SphSymbolic /. subsVal]]]];
Print["Cartesian-form check -- max |Cart2Sph(x,y,z) - Cart2Sph(theta,phi)| = ",
  cartMaxDiff, " (should be ~0)"];

(* ============================================================ *)
(* ADS equation string conversion                                *)
(* ============================================================ *)

AdsString[expr_] := Module[{s},
   s = ToString[expr, InputForm];
   s = StringReplace[s, {
      "Sqrt["    -> "sqrt(",
      "Abs["     -> "abs(",
      "Cos["     -> "cos(",
      "Sin["     -> "sin(",
      "Tan["     -> "tan(",
      "Csc["     -> "csc(",
      "Sec["     -> "sec(",
      "Cot["     -> "cot(",
      "ArcCos["  -> "acos(",
      "ArcSin["  -> "asin(",
      "ArcTan["  -> "atan(",
      "ArcCsc["  -> "acsc(",
      "ArcSec["  -> "asec(",
      "ArcCot["  -> "acot(",
      "Cosh["    -> "cosh(",
      "Sinh["    -> "sinh(",
      "Tanh["    -> "tanh(",
      "Csch["    -> "csch(",
      "Sech["    -> "sech(",
      "Coth["    -> "coth(",
      "ArcCosh[" -> "acosh(",
      "ArcSinh[" -> "asinh(",
      "ArcTanh[" -> "atanh(",
      "Pi" -> "pi",
      "[" -> "(",
      "]" -> ")",
      "^" -> "**"
      }];
   s = StringReplace[s, Whitespace -> ""];
   s
   ];

(* ============================================================ *)
(* Shared row-building helper: (i,j) in {1,2,3}x{1,2,3}           *)
(* ============================================================ *)

BuildMatrixIndexedEntries[computeFn_] :=
  Flatten[Table[{i, j, computeFn[i, j]}, {i, 1, 3}, {j, 1, 3}], 1];

(* ============================================================ *)
(* ADS equation-file builder                                     *)
(* ============================================================ *)

(* Cart2Sph_ij: angle form, function of theta,phi, used directly at the
   observation direction (a direct input, no inverse trig needed).

   M_ij, r_of_xyz, rho_of_xyz: generic Cartesian-ratio template in x,y,z,r,
   not evaluated at any point. The Python script uses M_ji with the
   body-frame variables x_prime,... to get Sph2Cart at (theta',phi')
   without acos/atan2 in the matrix entries. The transpose holds at the
   same point only: Cart2Sph at the observation point and Sph2Cart at the
   body-frame point are in general not transposes of each other.

   theta_of_xyz, phi_of_xyz (theta_prime, phi_prime in ADS) give the
   body-frame sampling direction for GradYlm_vars_prime; M_ij does not use
   them. theta_of_xyz = acos(z/r) is already in [0,pi]; phi_of_xyz =
   atan2(y,x) is wrapped from (-pi,pi] to [0,2*pi) with ADS
   "if ... then ... else ... endif". The two phi lines are written as
   literal ADS text because AdsString does not map two-argument ArcTan or
   produce the if syntax. *)
BuildADSEquationFile[] := Module[{cart2sphAngleRows, bodyLinesCart2SphAngle,
   genericRows, bodyLinesGeneric, rLine, rhoLine, thetaLine, phiRawLine, phiLine},
   (* ---- Cart2Sph_ij: plain angle-based trig, function of theta,phi ---- *)
   cart2sphAngleRows = BuildMatrixIndexedEntries[GetCart2Sphij];
   bodyLinesCart2SphAngle = Table[
     "Cart2Sph_" <> ToString[row[[1]]] <> ToString[row[[2]]] <> "=" <>
      AdsString[row[[3]]] <> ";",
     {row, cart2sphAngleRows}
     ];

   (* ---- M_ij,r_of_xyz,rho_of_xyz,theta_of_xyz,phi_of_xyz: GENERIC
      Cartesian-ratio template, function of plain x,y,z,r -- NOT
      evaluated at any specific point ---- *)
   genericRows = BuildMatrixIndexedEntries[GetCart2SphCartesianij];
   bodyLinesGeneric = Table[
     "M_" <> ToString[row[[1]]] <> ToString[row[[2]]] <> "=" <>
      AdsString[row[[3]]] <> ";",
     {row, genericRows}
     ];
   rLine = "r_of_xyz=" <> AdsString[rDef] <> ";";
   rhoLine = "rho_of_xyz=" <> AdsString[rhoDef] <> ";";
   thetaLine = "theta_of_xyz=" <> AdsString[ArcCos[z/r]] <> ";";
   phiRawLine = "phi_of_xyz_raw=atan2(y,x);";
   phiLine = "phi_of_xyz=if (phi_of_xyz_raw<0) then (phi_of_xyz_raw+2*pi) else (phi_of_xyz_raw) endif;";
   StringJoin[
    "// Cart2Sph_ij: PLAIN ANGLE-BASED trig, function of theta,phi --\n",
    "// use directly for the OBSERVATION direction (a direct input, no\n",
    "// inverse trig needed).\n",
    "//\n",
    "// M_ij,r_of_xyz,rho_of_xyz,theta_of_xyz,phi_of_xyz: GENERIC\n",
    "// Cartesian-ratio template, function of plain x,y,z,r, NOT\n",
    "// evaluated at any specific point. create_vector_field_rotation_vars.py\n",
    "// uses this, index-transposed (M_ij -> M_ji) and with x,y,z,r,rho\n",
    "// renamed to x_prime,y_prime,z_prime,r_prime,rho_prime, to build Sph2Cart at\n",
    "// the BODY-FRAME direction -- avoiding any inverse trig on\n",
    "// theta_prime,phi_prime inside the matrix entries themselves.\n",
    "// theta_of_xyz,phi_of_xyz (renamed theta_prime,phi_prime) are\n",
    "// REPORTING-ONLY there -- not used inside M_ij. theta_of_xyz in\n",
    "// [0,pi] via ArcCos's principal range (no correction needed);\n",
    "// phi_of_xyz wrapped into [0,2*pi) since atan2 alone returns\n",
    "// (-pi,pi].\n",
    "//\n",
    "// Paste into an ADS Eqn (equation) component; theta,phi (for\n",
    "// Cart2Sph_ij) and x,y,z (for the generic template) must already\n",
    "// be defined as ADS variables -- see create_vector_field_rotation_vars.py.\n",
    StringRiffle[bodyLinesCart2SphAngle, "\n"], "\n",
    rLine, "\n",
    rhoLine, "\n",
    StringRiffle[bodyLinesGeneric, "\n"], "\n",
    thetaLine, "\n",
    phiRawLine, "\n",
    phiLine
    ]
   ];

(* ============================================================ *)
(* In-notebook preview: nicely typeset Grid                       *)
(* ============================================================ *)

BuildPreviewGrid[computeFn_] := Module[{headerRow, bodyRows},
   headerRow = Prepend[{1, 2, 3}, "i/j"];
   bodyRows = Table[
     Prepend[Table[computeFn[i, j] // TraditionalForm, {j, 1, 3}], i],
     {i, 1, 3}
     ];
   Grid[
    Prepend[bodyRows, headerRow],
    Frame -> All,
    Background -> {None, {LightGray, None}},
    Alignment -> Center,
    Spacings -> {2, 1.5}
    ]
   ];

Print["Preview: Cart2Sph(theta,phi) -- angle form, used directly for the \
OBSERVATION direction"];
Print[BuildPreviewGrid[GetCart2Sphij]];
Print["Preview: M(x,y,z,r,rho) -- generic Cartesian-ratio template, used \
transposed+renamed for Sph2Cart at the BODY-FRAME direction"];
Print[BuildPreviewGrid[GetCart2SphCartesianij]];

(* ---- Build and export the single ADS equation file ---- *)
adsEquations = BuildADSEquationFile[];
adsFile = Export[FileNameJoin[{"tables", "Cart2Sph_ADS_Equations.txt"}], adsEquations, "Text"];
Print["Cart2Sph ADS equation file written to: ", adsFile];

Print["All exports complete."];
