(* rotation3D.wl
   Symbolic 3x3 rotation matrix R(alpha, Ux, Uy, Uz) (Rodrigues' formula:
   rotation by alpha [rad] about the axis (Ux,Uy,Uz)), checked against
   RotationMatrix[]. Writes to tables/:
     Rotation_table.txt          tab-delimited MATLAB str2sym strings
     Rotation_table.tex          LaTeX pmatrix
     Rotation_ADS_Equations.txt  ADS Eqn text (magU, R11..R33; alpha -> rot_angle) *)

ClearAll[magUDef, RotationMatrixSymbolic, GetRij, MatlabString, AdsString,
  LatexEntryTall, BuildRotationIndexedEntries, BuildRotationMatrixLatex,
  BuildADSEquationFile, BuildPreviewGrid];

(* Run from the mathematica/ folder; generated tables are written to tables/. *)
SetDirectory[NotebookDirectory[]];
If[!DirectoryQ["tables"], CreateDirectory["tables"]];
Print["Working directory set to: ", Directory[]];

(* ============================================================ *)
(* Core symbolic definition                                      *)
(* ============================================================ *)

(* Rotation matrix R(alpha,Ux,Uy,Uz) by Rodrigues' formula: R.x rotates x
   by alpha radians about the line through the origin and (Ux,Uy,Uz). The
   axis is normalised internally, element by element as in rotationmat3D.m
   (Bileschi 2009). alpha, Ux, Uy, Uz are free symbols with the same names
   as the ADS parameters.

   The axis magnitude magU = sqrt(Ux^2+Uy^2+Uz^2) is kept as a free symbol
   (u = Ux/magU, ...), so each entry simplifies to a rational function of
   Ux,Uy,Uz,magU and trig of alpha. magUDef holds its definition, used for
   the numerical check and written as a separate equation in the exports. *)
magUDef = Sqrt[Ux^2 + Uy^2 + Uz^2];

(* Module keeps the helper symbols local. *)
RotationMatrixSymbolic = Module[{u, v, w, u2, v2, w2, c, s},
   u = Ux/magU;
   v = Uy/magU;
   w = Uz/magU;
   u2 = u^2;
   v2 = v^2;
   w2 = w^2;
   c = Cos[alpha];
   s = Sin[alpha];
   Simplify[
    {
     {u2 + (v2 + w2)*c, u*v*(1 - c) - w*s, u*w*(1 - c) + v*s},
     {u*v*(1 - c) + w*s, v2 + (u2 + w2)*c, v*w*(1 - c) - u*s},
     {u*w*(1 - c) - v*s, v*w*(1 - c) + u*s, w2 + (u2 + v2)*c}
     },
    Assumptions -> {alpha \[Element] Reals, Ux \[Element] Reals,
      Uy \[Element] Reals, Uz \[Element] Reals, magU > 0}
    ]
   ];

(* Wrapper: returns the (i,j) entry symbolically, i,j in {1,2,3} *)
GetRij[i_, j_] := RotationMatrixSymbolic[[i, j]];

(* ============================================================ *)
(* Numerical cross-check against Mathematica's built-in           *)
(* RotationMatrix[] (independent of rotationmat3D.m) -- should     *)
(* agree to machine precision.                                     *)
(* ============================================================ *)

SeedRandom[42];
alphaVal = RandomReal[{-Pi, Pi}];
UxVal = RandomReal[{-2, 2}];
UyVal = RandomReal[{-2, 2}];
UzVal = RandomReal[{-2, 2}];
(* magU must be evaluated from magUDef at the same Ux,Uy,Uz. *)
magUVal = magUDef /. {Ux -> UxVal, Uy -> UyVal, Uz -> UzVal};
testVals = {
   alpha -> alphaVal, Ux -> UxVal, Uy -> UyVal, Uz -> UzVal,
   magU -> magUVal
   };
numericOurs = N[RotationMatrixSymbolic /. testVals];
numericBuiltin = N[RotationMatrix[alpha, {Ux, Uy, Uz}] /. testVals];
maxDiff = Max[Abs[Flatten[numericOurs - numericBuiltin]]];
Print["Sanity check -- max |our R - built-in RotationMatrix[]| = ",
  maxDiff, " (should be ~0)"];

(* ============================================================ *)
(* MATLAB string conversion                                      *)
(* ============================================================ *)

(* MATLAB/str2sym-compatible string, wrapped in double quotes so MATLAB's
   text import reads it as one field (same mapping as LegendreFThetaFPhi.wl). *)
MatlabString[expr_] := Module[{s},
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
      "]" -> ")"
      }];
   "\"" <> s <> "\""
   ];

(* ============================================================ *)
(* LaTeX table conversion                                        *)
(* ============================================================ *)

(* Plain TeXForm of each entry, for the LaTeX export only: alpha ->
   \[Alpha] on the expression; Ux,Uy,Uz,magU -> u_{x},u_{y},u_{z},
   u_{\mathrm{mag}} on the TeXForm string, which writes multi-letter
   symbols either bare or as \text{...}. *)
LatexEntryTall[expr_] := Module[{texExpr, s},
   texExpr = expr /. {alpha -> \[Alpha]};
   s = ToString[TeXForm[texExpr]];
   s = StringReplace[s, {
      "\\text{Ux}"   -> "u_{x}",
      "\\text{Uy}"   -> "u_{y}",
      "\\text{Uz}"   -> "u_{z}",
      "\\text{magU}" -> "u_{\\mathrm{mag}}",
      "Ux"   -> "u_{x}",
      "Uy"   -> "u_{y}",
      "Uz"   -> "u_{z}",
      "magU" -> "u_{\\mathrm{mag}}"
      }];
   s
   ];

(* ============================================================ *)
(* ADS equation string conversion                                *)
(* ============================================================ *)

(* ADS-equation string: same mapping as MatlabString, unquoted, without
   whitespace, and ^ written as ** (ADS power operator). *)
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

BuildRotationIndexedEntries[computeFn_] :=
  Flatten[Table[{i, j, computeFn[i, j]}, {i, 1, 3}, {j, 1, 3}], 1];

(* ============================================================ *)
(* LaTeX builder: the 9 entries as a 3x3 pmatrix (rows i, cols j) *)
(* ============================================================ *)

BuildRotationMatrixLatex[computeFn_, captionLabel_, entryFn_] :=
  Module[{rowStrs},
   rowStrs = Table[
     StringRiffle[Table[entryFn[computeFn[i, j]], {j, 1, 3}], " & "],
     {i, 1, 3}
     ];
   StringJoin[
    "\\begin{center}\n",
    "\\textbf{$", captionLabel, "$}\\\\[6pt]\n",
    "\\[\n",
    "R(\\alpha, u_{x}, u_{y}, u_{z}) =\n",
    "\\begin{pmatrix}\n",
    StringRiffle[rowStrs, " \\\\\n"], "\n",
    "\\end{pmatrix}\n",
    "\\]\n",
    "where $u_{\\mathrm{mag}} = \\sqrt{u_x^{2}+u_y^{2}+u_z^{2}}$.\n",
    "\\end{center}"
    ]
   ];

(* ============================================================ *)
(* ADS equation-file builder                                     *)
(* ============================================================ *)

(* ADS equations for an Eqn component (R depends on continuous parameters,
   so no Discrete lookup file): magU first, then R11..R33. rot_angle, Ux,
   Uy, Uz must be defined upstream (e.g. in a VAR block). alpha is renamed
   rot_angle in this export only, to avoid clashing with other ADS names. *)
BuildADSEquationFile[computeFn_, matName_] := Module[{rows, magULine, bodyLines},
   rows = BuildRotationIndexedEntries[computeFn];
   magULine = "magU=" <> AdsString[magUDef] <> ";";
   bodyLines = Table[
     matName <> ToString[row[[1]]] <> ToString[row[[2]]] <> "=" <>
      AdsString[row[[3]] /. {alpha -> rot_angle}] <> ";",
     {row, rows}
     ];
   StringJoin[
    "// Rotation matrix R(rot_angle,Ux,Uy,Uz), Rodrigues' rotation formula\n",
    "// magU = |axis vector| = sqrt(Ux^2+Uy^2+Uz^2), computed once here\n",
    "// and reused below rather than recomputed in every Rij entry.\n",
    "// Paste into an ADS Eqn (equation) component; rot_angle, Ux, Uy, Uz\n",
    "// must already be defined as ADS variables (e.g. in a VAR block).\n",
    magULine, "\n",
    StringRiffle[bodyLines, "\n"]
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

Print["Preview: R(alpha,Ux,Uy,Uz,magU)"];
Print[BuildPreviewGrid[GetRij]];

(* ---- Build MATLAB-string table (3x3, rows=i, columns=j) ---- *)
(* First row: magU definition (entries reference magU); columns 3-4 are
   empty padding. *)
RotationTable = {
   {"magU_def", MatlabString[magUDef], "", ""},
   Prepend[{1, 2, 3}, "i/j"]
   };
Do[
  Module[{rowR},
   rowR = {iVal};
   Do[AppendTo[rowR, MatlabString[GetRij[iVal, jVal]]], {jVal, 1, 3}];
   AppendTo[RotationTable, rowR];
   ],
  {iVal, 1, 3}
  ];

(* ---- Export tab-delimited text file (for MATLAB str2sym import) ---- *)
RotationText = StringRiffle[
   Map[StringRiffle[ToString /@ #, "\t"] &, RotationTable],
   "\n"
   ];
txtFile = Export[FileNameJoin[{"tables", "Rotation_table.txt"}], RotationText, "Text"];
Print["Rotation matrix TXT (tab-delimited) written to: ", txtFile];

(* ---- Build and export LaTeX table ---- *)
latexRotation = BuildRotationMatrixLatex[GetRij, "R(\\alpha,u_{x},u_{y},u_{z})",
   LatexEntryTall];
texFile = Export[FileNameJoin[{"tables", "Rotation_table.tex"}], latexRotation, "Text"];
Print["Rotation matrix LaTeX table written to: ", texFile];

(* ---- Build and export ADS equation file ---- *)
adsRotation = BuildADSEquationFile[GetRij, "R"];
adsFile = Export[FileNameJoin[{"tables", "Rotation_ADS_Equations.txt"}], adsRotation, "Text"];
Print["Rotation matrix ADS equation file written to: ", adsFile];

Print["All exports complete."];
