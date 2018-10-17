(* ::Package:: *)

(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)

(* :Title: longskindepthewjn *)
(* :Context: longskindepthewjn` *)
(* :Author: Deepak *)
(* :Date: 2018-08-15 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: *)
(* :Copyright: (c) 2018 Deepak *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["longskindepthewjn`"];
(* Exported symbols added here with SymbolName::usage *)

Unprotect["longskindepthewjn`*"];
ClearAll["longskindepthewjn`*"];
ClearAll["longskindepthewjn`Private`*"];

unscaledJ::usage = "unscaledJ[aSource, region, zeroPoint, nHat] returns J[x,y,z] for some 3D region region without prefactors";
unscaledBz::usage = "unscaledB[j, radius] returns the magnetic field without prefactors for a current j defined over cube bounded by radius r";


Begin["`Private`"];

solveForF3D[aSource_, region_, zeroPointPred_, nHat_] := NDSolveValue[{
	Laplacian[f[x, y, z], {x, y, z}] == NeumannValue[aSource[x, y, z] . nHat[x, y, z], True],
	DirichletCondition[f[x, y, z] == 0, zeroPointPred[x, y, z]]
}, f, Element[{x, y, z}, region]];

(* Utility to simplify some of the type handling here, to make sure we can move the piecewise inside the list for j*)
piecewise3DFunc[function_, region_, default_, index_] := Function[{x, y, z},
		Simplify`PWToUnitStep[Piecewise[
			{{function[x, y, z][[index]], {x, y, z} \[Element] region}}
			, default
		]]
	];

unscaledJ[aSource_, region_, zeroPoint_, nHat_] := Module[{f, gradFPlusAs, j},
	f = solveForF3D[aSource, region, zeroPoint, nHat];
	gradFPlusAs[x_, y_, z_] := Evaluate[D[f[x, y, z], {{x, y, z}}] + aSource[x, y, z]];
	j[x_, y_, z_] = Through[
			Thread[Unevaluated[piecewise3DFunc[gradFPlusAs, region, {0, 0, 0}, {1, 2, 3}]]]
			[x, y, z]
		];
	j
];

unscaledBz[j_, radius_] := Module[{jSafe, bz},
	jSafe[x_, y_, z_] := Quiet[Check[j[x, y, z], {0, 0, 0}, InterpolatingFunction::femdmval], InterpolatingFunction::femdmval];
	Print[jSafe[1, 0, 0]];
	Print[jSafe[10, 0, 0]];
	Print[TemplateApply[
		StringTemplate["Creating bz with radius `1`"],
		radius
	]];
	bz[x_, y_, z_] := NIntegrate[
		((jSafe[xp, yp, zp][[1]]*(y - yp)) - (jSafe[xp, yp, zp][[2]]*(x - xp)) )/(Norm[{x, y, z} - {xp, yp, zp}]^3),
		{xp, -radius, radius}, {yp, -radius, radius}, {zp, -radius, radius},
		WorkingPrecision -> 100
	];
	bz
];

End[]; (* `Private` *)
Protect["longskindepthewjn`*"];
EndPackage[];
