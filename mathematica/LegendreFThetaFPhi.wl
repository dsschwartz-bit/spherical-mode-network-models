(* ::Package:: *)

(* LegendreFThetaFPhi.wl
   Symbolic angular functions of the spherical modes, built from the fully
   normalised associated Legendre function Pbar_lm(x), x = cos(th):
     F_theta(l,m) = -sqrt(1-x^2) dPbar_lm/dx   (= -sin(th) dPbar/dx)
     F_phi(l,m)   =  Pbar_lm(x) / sqrt(1-x^2)  (= Pbar/sin(th))
   for l = 0..lMax, m = -l..l. Writes to tables/:
     FTheta_table.txt, FPhi_table.txt                  tab-delimited l x m tables of MATLAB
                                                       str2sym strings (variable th)
     FTheta_table.tex, FPhi_table.tex                  LaTeX longtables
     FTheta_ADS_Discrete.txt, FPhi_ADS_Discrete.txt    ADS DataAccessComponent "Discrete"
                                                       files, index i = l(l+1)+m (variable theta) *)

ClearAll[NormFactor, NormAssocLegendre, NormAssocLegendreDeriv,
  ComputeFThetaFPhi, MatlabString, LatexEntry, BuildLatexTable,
  GetFTheta, GetFPhi];

(* Run from the mathematica/ folder; generated tables are written to tables/. *)
SetDirectory[NotebookDirectory[]];
If[!DirectoryQ["tables"], CreateDirectory["tables"]];
Print["Working directory set to: ", Directory[]];

(* ============================================================ *)
(* Core symbolic definitions                                     *)
(* ============================================================ *)

(* Normalisation factor without the 2*Pi of the phi integration, including
   (-1)^m to cancel the Condon-Shortley phase of LegendreP *)
NormFactor[l_, m_] := (-1)^m * Sqrt[(2 l + 1)/2 * (l - m)!/(l + m)!];

(* Fully normalized associated Legendre function *)
NormAssocLegendre[l_, m_, x_] := NormFactor[l, m] * LegendreP[l, m, x];

(* Derivative of the normalized function w.r.t. x *)
NormAssocLegendreDeriv[l_, m_, x_] := D[NormAssocLegendre[l, m, x], x];

(* Main function: computes F_theta and F_phi for given l, m
   (symbolic, in terms of th = theta) *)
ComputeFThetaFPhi[l_, m_] := Module[
   {exprP, exprDP, FThetaExpr, FPhiExpr, FTheta, FPhi},
   exprP = NormAssocLegendre[l, m, x];
   exprDP = NormAssocLegendreDeriv[l, m, x];

   (* F_theta: -Sqrt[1-x^2] * dP/dx, simplified, then x -> Cos[th] *)
   FThetaExpr = -Sqrt[1 - x^2] * exprDP // Together // Cancel // Simplify;
   FTheta = FThetaExpr /. x -> Cos[th] //
      Simplify[#, 0 <= th <= Pi] & // TrigReduce // Simplify;

   (* F_phi: (1/Sqrt[1-x^2]) * P, simplified, then x -> Cos[th] *)
   FPhiExpr = (1/Sqrt[1 - x^2]) * exprP // Together // Cancel // Simplify;
   FPhi = FPhiExpr /. x -> Cos[th] //
      Simplify[#, 0 <= th <= Pi] & // TrigReduce // Simplify;

   {FTheta, FPhi}
   ];

(* Wrapper functions returning just F_theta or F_phi symbolically *)
GetFTheta[l_, m_] := First[ComputeFThetaFPhi[l, m]];
GetFPhi[l_, m_] := Last[ComputeFThetaFPhi[l, m]];

(* ============================================================ *)
(* MATLAB string conversion                                      *)
(* ============================================================ *)

(* MATLAB/str2sym-compatible string of an expression, wrapped in double
   quotes so MATLAB's text import reads it as a single field. *)
MatlabString[expr_] := Module[{s},
   s = ToString[expr, InputForm];
   s = StringReplace[s, {
      "Sqrt["    -> "sqrt(",
      "Abs["     -> "abs(",
      (* Basic trig *)
      "Cos["     -> "cos(",
      "Sin["     -> "sin(",
      "Tan["     -> "tan(",
      "Csc["     -> "csc(",
      "Sec["     -> "sec(",
      "Cot["     -> "cot(",
      (* Inverse trig *)
      "ArcCos["  -> "acos(",
      "ArcSin["  -> "asin(",
      "ArcTan["  -> "atan(",
      "ArcCsc["  -> "acsc(",
      "ArcSec["  -> "asec(",
      "ArcCot["  -> "acot(",
      (* Hyperbolic *)
      "Cosh["    -> "cosh(",
      "Sinh["    -> "sinh(",
      "Tanh["    -> "tanh(",
      "Csch["    -> "csch(",
      "Sech["    -> "sech(",
      "Coth["    -> "coth(",
      (* Inverse hyperbolic *)
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

(* The LaTeX helpers replace th/ph by \[Theta]/\[Phi] before TeXForm so the
   tables show \theta/\phi; the MATLAB and ADS exports keep th/theta. *)

(* Plain TeXForm (fractions left as TeXForm writes them). Alternative
   entry formatter; both tables currently use LatexEntryCoeff. *)
LatexEntryTall[expr_] := Module[{texExpr},
   texExpr = expr /. {th -> \[Theta], ph -> \[Phi]};
   ToString[TeXForm[texExpr]]
   ];

(* Writes the rational factor of a product (e.g. -1/4) as a compact
   \frac{}{} in front of the remaining factors, instead of one tall
   fraction. Irrational factors such as Sqrt[15] stay in the rest. *)
LatexEntryCoeff[expr_] := Module[
   {texExpr, factors, ratFactors, restFactors, ratCoeff, restExpr,
    sign, absCoeff, numStr, denStr, restStr},
   texExpr = expr /. {th -> \[Theta], ph -> \[Phi]};
   factors = If[Head[texExpr] === Times, List @@ texExpr, {texExpr}];
   ratFactors = Select[factors, (IntegerQ[#] || Head[#] === Rational) &];
   restFactors = Complement[factors, ratFactors];
   ratCoeff = If[ratFactors === {}, 1, Times @@ ratFactors];
   restExpr = If[restFactors === {}, 1, Times @@ restFactors];
   Which[
    restExpr === 1, ToString[TeXForm[ratCoeff]],
    ratCoeff === 1, ToString[TeXForm[restExpr]],
    True,
      sign = If[ratCoeff < 0, "-", ""];
      absCoeff = Abs[ratCoeff];
      numStr = ToString[Numerator[absCoeff]];
      denStr = ToString[Denominator[absCoeff]];
      (* Bracket a sum so the coefficient multiplies all of it. *)
      restStr = If[Head[restExpr] === Plus,
        "\\left(" <> ToString[TeXForm[restExpr]] <> "\\right)",
        ToString[TeXForm[restExpr]]
        ];
      sign <> "\\frac{" <> numStr <> "}{" <> denStr <> "} " <> restStr
    ]
   ];

(* Rows {i, l, m, value} sorted by the linear index i = l(l+1)+m (i = 0
   for l = m = 0). Used for the LaTeX row order and the ADS index column. *)
BuildIndexedRows[lValues_, computeFn_] := SortBy[
   Flatten[
    Table[
     If[Abs[mVal] <= lVal,
      {lVal*(lVal + 1) + mVal, lVal, mVal, computeFn[lVal, mVal]},
      Nothing
      ],
     {lVal, lValues}, {mVal, -lVal, lVal}
     ],
    1],
   First
   ];

(* Two-column longtable: (l,m) label and value, one row per (l,m), so the
   width does not grow with lMax and the table breaks across pages. *)
BuildLatexTable[lValues_, lMax_, computeFn_, captionLabel_, headerLabel_, entryFn_] := Module[
   {allRows, header, bodyRows},
   allRows = BuildIndexedRows[lValues, computeFn];
   header = "\\boldmath \\rule[-9pt]{0pt}{21pt}$(l,m)$ & " <> headerLabel <> "\\unboldmath \\\\";
   bodyRows = Table[
     Module[{row = allRows[[i]], rowColorStr, entryStr},
      (* \rowcolor per row: \rowcolors is unreliable with longtable. *)
      rowColorStr = If[EvenQ[i], "\\rowcolor{gray!15}\n", ""];
      entryStr = StringRiffle[{
         "$(l=" <> ToString[row[[2]]] <> ", m=" <> ToString[row[[3]]] <> ")$",
         "$" <> entryFn[row[[4]]] <> "$"
         }, " & "] <> " \\\\ \\hline";
      rowColorStr <> entryStr
      ],
     {i, Length[allRows]}
     ];
   StringJoin[
    "\\begin{center}\n",
    "\\begin{longtable}{|c|c|}\n",
    "\\caption{", captionLabel, "} \\\\\n",
    "\\hline\n", header, "\n\\hline\n",
    "\\endfirsthead\n",
    "\\hline\n", header, "\n\\hline\n",
    "\\endhead\n",
    StringRiffle[bodyRows, "\n"], "\n",
    "\\end{longtable}\n",
    "\\end{center}"
    ]
   ];

(* ============================================================ *)
(* Keysight ADS DataAccessComponent "Discrete" file export       *)
(* ============================================================ *)

(* ADS-equation string of an expression. Same function mapping as
   MatlabString, but unquoted (Discrete files hold string columns
   unquoted), with all whitespace removed (Discrete columns are
   whitespace-delimited) and ^ written as ** (ADS power operator). *)
AdsString[expr_] := Module[{texExpr, s},
   (* Rename th -> theta on the expression, not the string, so names
      such as coth are not altered. *)
   texExpr = expr /. {th -> theta};
   s = ToString[texExpr, InputForm];
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

(* ADS "Discrete" data file (BEGIN DSCRDATA ... END DSCRDATA): index
   column i = l(l+1)+m in ascending order (matches ADS "Index Lookup",
   0-based row position) and the formula string as dependent variable. *)
BuildADSDiscreteFile[lValues_, computeFn_, depVarName_] := Module[
   {rows, bodyLines},
   rows = BuildIndexedRows[lValues, computeFn];
   bodyLines = Table[
     ToString[row[[1]]] <> " " <> AdsString[row[[4]]],
     {row, rows}
     ];
   StringJoin[
    (* REM line recording lMax, read by the Python scripts. *)
    "REM lMax=", ToString[Max[lValues]], "\n",
    "BEGIN DSCRDATA\n",
    "% INDEX ", depVarName, "\n",
    StringRiffle[bodyLines, "\n"], "\n",
    "END DSCRDATA"
    ]
   ];

(* ============================================================ *)
(* Iterate over given l values; m ranges from -l to l            *)
(* ============================================================ *)

lValues = Range[0, 11];
lMax = Max[lValues];

(* ============================================================ *)
(* In-notebook preview: typeset l x m grid of the expressions     *)
(* ============================================================ *)

BuildPreviewGrid[lValues_, lMax_, computeFn_] := Module[
   {headerRow, bodyRows},
   headerRow = Prepend[
     Table[m, {m, -lMax, lMax}],
     "l/m"
     ];
   bodyRows = Table[
     Prepend[
      Table[
       If[Abs[mVal] <= lVal, computeFn[lVal, mVal] // TraditionalForm, 0],
       {mVal, -lMax, lMax}
       ],
      lVal
      ],
     {lVal, lValues}
     ];
   Grid[
    Prepend[bodyRows, headerRow],
    Frame -> All,
    Background -> {None, {LightGray, None}},
    Alignment -> Center,
    Spacings -> {2, 1.5}
    ]
   ];

Print["Preview: F_\[Theta](l,m)"];
Print[BuildPreviewGrid[lValues, lMax, GetFTheta]];

Print["Preview: F_\[Phi](l,m)"];
Print[BuildPreviewGrid[lValues, lMax, GetFPhi]];

(* ---- MATLAB tables: one row per l, columns m = -lMax..lMax; first cell "l/m" ---- *)
mHeaderCols = Table[mVal, {mVal, -lMax, lMax}];
FThetaTable = {Prepend[mHeaderCols, "l/m"]};
FPhiTable = {Prepend[mHeaderCols, "l/m"]};

Do[
  Module[{rowFTheta, rowFPhi, FThetaResult, FPhiResult},
   rowFTheta = {lVal};
   rowFPhi = {lVal};
   Do[
     If[Abs[mVal] <= lVal,
      {FThetaResult, FPhiResult} = ComputeFThetaFPhi[lVal, mVal];
      AppendTo[rowFTheta, MatlabString[FThetaResult]];
      AppendTo[rowFPhi, MatlabString[FPhiResult]];
      ,
      (* m out of range for this l: entry is identically 0 *)
      AppendTo[rowFTheta, MatlabString[0]];
      AppendTo[rowFPhi, MatlabString[0]];
      ],
     {mVal, -lMax, lMax}
     ];
   AppendTo[FThetaTable, rowFTheta];
   AppendTo[FPhiTable, rowFPhi];
   ],
  {lVal, lValues}
  ];

(* ---- Export tab-delimited text files (for MATLAB str2sym import) ----
   Joined manually and exported as "Text" (not "TSV") so the already
   quoted expression fields are not escaped again. *)
FThetaText = StringRiffle[
   Map[StringRiffle[ToString /@ #, "\t"] &, FThetaTable],
   "\n"
   ];
FPhiText = StringRiffle[
   Map[StringRiffle[ToString /@ #, "\t"] &, FPhiTable],
   "\n"
   ];

txtFile1 = Export[FileNameJoin[{"tables", "FTheta_table.txt"}], FThetaText, "Text"];
txtFile2 = Export[FileNameJoin[{"tables", "FPhi_table.txt"}], FPhiText, "Text"];

Print["F_theta TXT (tab-delimited) written to: ", txtFile1];
Print["F_phi TXT (tab-delimited) written to: ", txtFile2];

(* ---- Build and export LaTeX tables ---- *)
latexFTheta = BuildLatexTable[lValues, lMax, GetFTheta, "$F_\\theta(l,m)$",
   "$F_\\theta := -\\sin\\theta \\, \\frac{d\\bar{P}_{l,m}}{dx}\\Big|_{x=\\cos\\theta}$",
   LatexEntryCoeff];
latexFPhi = BuildLatexTable[lValues, lMax, GetFPhi, "$F_\\phi(l,m)$",
   "$F_\\phi := \\frac{\\bar{P}_{l,m}(\\cos\\theta)}{\\sin\\theta}$",
   LatexEntryCoeff];

texFile1 = Export[FileNameJoin[{"tables", "FTheta_table.tex"}], latexFTheta, "Text"];
texFile2 = Export[FileNameJoin[{"tables", "FPhi_table.tex"}], latexFPhi, "Text"];

Print["F_theta LaTeX table written to: ", texFile1];
Print["F_phi LaTeX table written to: ", texFile2];

(* ---- Build and export ADS DAC "Discrete" data files ---- *)
adsFTheta = BuildADSDiscreteFile[lValues, GetFTheta, "Ftheta"];
adsFPhi = BuildADSDiscreteFile[lValues, GetFPhi, "Fphi"];

adsFile1 = Export[FileNameJoin[{"tables", "FTheta_ADS_Discrete.txt"}], adsFTheta, "Text"];
adsFile2 = Export[FileNameJoin[{"tables", "FPhi_ADS_Discrete.txt"}], adsFPhi, "Text"];

Print["F_theta ADS Discrete file written to: ", adsFile1];
Print["F_phi ADS Discrete file written to: ", adsFile2];

Print["All exports complete."];
