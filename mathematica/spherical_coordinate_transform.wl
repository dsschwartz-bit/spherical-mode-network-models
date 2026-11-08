(* spherical_coordinate_transform.wl
   Spherical <-> Cartesian vector-component matrices for rotating a field
   sample together with its sample direction (rotated position
   (x_rot,y_rot,z_rot) = R.(x,y,z)). Writes to tables/:
     Sph2Cart_Cart2Sph_ADS_Equations.txt  ADS Eqn text: Sph2Cart_ij(theta,phi),
                                          Cart2Sph_ij(x_rot,y_rot,z_rot),
                                          theta_rot, phi_rot
   The python/ads scripts use the observation-point formulation of
   vector_field_spherical_coordinate_transform.wl instead. *)

ClearAll[Cart2SphSymbolic, Sph2CartSymbolic, Cart2SphCartesianSymbolic,
  rDef, rhoDef, GetCart2Sphij, GetSph2Cartij, GetCart2SphCartesianij,
  AdsString, BuildMatrixIndexedEntries, BuildSphCartADSEquationFile,
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

(* Cart2Sph is needed at the rotated sample direction (x_rot,y_rot,z_rot)
   = R.(x,y,z), with (x,y,z) = rhat(theta,phi) rotated by the same matrix
   as the field vector. To avoid inverse trig, the trig functions are
   written as ratios of x,y,z (r = sqrt(x^2+y^2+z^2), rho = sqrt(x^2+y^2)):

     sin(theta) = rho/r     cos(theta) = z/r
     sin(phi)   = y/rho     cos(phi)   = x/rho

   Substituting these into Cart2SphSymbolic and simplifying gives:

     rhat     = ( x/r,          y/r,          z/r   )
     thetahat = ( x*z/(r*rho),  y*z/(r*rho), -rho/r  )
     phihat   = (-y/rho,        x/rho,        0      )

   r and rho are kept as free symbols (definitions rDef, rhoDef), as magU
   in rotation3D.wl, so the ADS export computes each square root once. *)
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

(* Sph2Cart_ij: function of theta,phi (the original sample direction).
   Cart2Sph_ij: Cartesian form with x,y,z,r,rho -> x_rot,y_rot,z_rot,
   r_rot,rho_rot; x_rot,y_rot,z_rot must be defined upstream in ADS (the
   position r*rhat(theta,phi) rotated by Rotation_lib:RotationMatrix_vars).
   No inverse trig appears in Cart2Sph_ij.

   theta_rot, phi_rot give the rotated sample location in spherical
   coordinates (outputs only). theta_rot = acos(z_rot/r_rot) is already in
   [0,pi]; phi_rot = atan2(y_rot,x_rot) is wrapped from (-pi,pi] to
   [0,2*pi) with ADS "if ... then ... else ... endif". The two phi_rot
   lines are written as literal ADS text because AdsString does not map
   two-argument ArcTan or produce the if syntax. *)
BuildSphCartADSEquationFile[] := Module[{sph2cartRows, cart2sphRows,
   bodyLinesS2C, rRotLine, rhoRotLine, bodyLinesC2S, thetaRotLine,
   phiRotRawLine, phiRotLine, rotSubs},
   rotSubs = {x -> x_rot, y -> y_rot, z -> z_rot, r -> r_rot, rho -> rho_rot};
   sph2cartRows = BuildMatrixIndexedEntries[GetSph2Cartij];
   cart2sphRows = BuildMatrixIndexedEntries[GetCart2SphCartesianij];
   bodyLinesS2C = Table[
     "Sph2Cart_" <> ToString[row[[1]]] <> ToString[row[[2]]] <> "=" <>
      AdsString[row[[3]]] <> ";",
     {row, sph2cartRows}
     ];
   rRotLine = "r_rot=" <> AdsString[rDef /. {x -> x_rot, y -> y_rot, z -> z_rot}] <> ";";
   rhoRotLine = "rho_rot=" <> AdsString[rhoDef /. {x -> x_rot, y -> y_rot, z -> z_rot}] <> ";";
   bodyLinesC2S = Table[
     "Cart2Sph_" <> ToString[row[[1]]] <> ToString[row[[2]]] <> "=" <>
      AdsString[row[[3]] /. rotSubs] <> ";",
     {row, cart2sphRows}
     ];
   thetaRotLine = "theta_rot=" <> AdsString[ArcCos[z_rot/r_rot]] <> ";";
   phiRotRawLine = "phi_rot_raw=atan2(y_rot,x_rot);";
   phiRotLine = "phi_rot=if (phi_rot_raw<0) then (phi_rot_raw+2*pi) else (phi_rot_raw) endif;";
   StringJoin[
    "// Sph2Cart: maps spherical vector components {Vr,Vtheta,Vphi},\n",
    "// specified at the ORIGINAL frame's (theta,phi), to Cartesian\n",
    "// components {Vx,Vy,Vz}.\n",
    "//\n",
    "// Cart2Sph: maps Cartesian components {Vx,Vy,Vz} back to spherical\n",
    "// components {Vr,Vtheta,Vphi}, evaluated at the ROTATED SAMPLE\n",
    "// LOCATION (x_rot,y_rot,z_rot) -- built purely from Cartesian\n",
    "// ratios (r_rot=sqrt(x_rot^2+y_rot^2+z_rot^2),\n",
    "// rho_rot=sqrt(x_rot^2+y_rot^2)), with NO sin/cos/inverse-trig of\n",
    "// theta_rot,phi_rot appearing in Cart2Sph_ij itself. x_rot,y_rot,\n",
    "// z_rot must already be defined upstream as the rotated position\n",
    "// (x,y,z)=r*rhat(theta,phi) after applying\n",
    "// Rotation_lib:RotationMatrix_vars -- this is a DIFFERENT rotated\n",
    "// vector than the field Vx_rot,Vy_rot,Vz_rot.\n",
    "//\n",
    "// theta_rot,phi_rot below are REPORTING-ONLY: the new sample\n",
    "// location's spherical coordinates, for output purposes -- not\n",
    "// used anywhere in Cart2Sph_ij above. theta_rot in [0,pi] via\n",
    "// ArcCos's principal range (no correction needed); phi_rot wrapped\n",
    "// into [0,2*pi) since atan2 alone returns (-pi,pi].\n",
    "//\n",
    "// Paste into an ADS Eqn (equation) component; theta, phi, x_rot,\n",
    "// y_rot, z_rot must already be defined as ADS variables.\n",
    StringRiffle[bodyLinesS2C, "\n"], "\n",
    rRotLine, "\n",
    rhoRotLine, "\n",
    StringRiffle[bodyLinesC2S, "\n"], "\n",
    thetaRotLine, "\n",
    phiRotRawLine, "\n",
    phiRotLine
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

Print["Preview: Sph2Cart(theta,phi)"];
Print[BuildPreviewGrid[GetSph2Cartij]];
Print["Preview: Cart2Sph(x,y,z,r,rho) -- the form actually exported for ADS \
(with x,y,z,r,rho -> x_rot,y_rot,z_rot,r_rot,rho_rot)"];
Print[BuildPreviewGrid[GetCart2SphCartesianij]];

(* ---- Build and export the single ADS equation file ---- *)
adsSphCart = BuildSphCartADSEquationFile[];
adsFile = Export[FileNameJoin[{"tables", "Sph2Cart_Cart2Sph_ADS_Equations.txt"}], adsSphCart, "Text"];
Print["Sph2Cart/Cart2Sph ADS equation file written to: ", adsFile];

Print["All exports complete."];
