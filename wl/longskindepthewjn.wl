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
unscaledBz::usage = "unscaledBz[j, radius] returns a function returning the magnetic field without prefactors in the z direction for a current j defined over cube bounded by radius r";
unscaledBy::usage = "unscaledBy[j, radius] returns a function returning the magnetic field without prefactors in the y direction for a current j defined over cube bounded by radius r";
unscaledBx::usage = "unscaledBx[j, radius] returns a function returning the magnetic field without prefactors in the x direction for a current j defined over cube bounded by radius r";


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
	bz[x_, y_, z_] := NIntegrate[
		((jSafe[xp, yp, zp][[1]]*(y - yp)) - (jSafe[xp, yp, zp][[2]]*(x - xp)) )/(Norm[{x, y, z} - {xp, yp, zp}]^3),
		{xp, -radius, radius}, {yp, -radius, radius}, {zp, -radius, radius},
		WorkingPrecision -> 100
	];
	bz
];

unscaledBx[j_, radius_] := Module[{jSafe, bx},
	jSafe[x_, y_, z_] := Quiet[Check[j[x, y, z], {0, 0, 0}, InterpolatingFunction::femdmval], InterpolatingFunction::femdmval];
	bx[x_, y_, z_] := NIntegrate[
		((jSafe[xp, yp, zp][[2]]*(z - zp)) - (jSafe[xp, yp, zp][[3]]*(y - yp)) )/(Norm[{x, y, z} - {xp, yp, zp}]^3),
		{xp, -radius, radius}, {yp, -radius, radius}, {zp, -radius, radius},
		WorkingPrecision -> 100
	];
	bx
];

unscaledBy[j_, radius_] := Module[{jSafe, by},
	jSafe[x_, y_, z_] := Quiet[Check[j[x, y, z], {0, 0, 0}, InterpolatingFunction::femdmval], InterpolatingFunction::femdmval];
	by[x_, y_, z_] := NIntegrate[
		((jSafe[xp, yp, zp][[3]]*(x - xp)) - (jSafe[xp, yp, zp][[1]]*(z - zp)) )/(Norm[{x, y, z} - {xp, yp, zp}]^3),
		{xp, -radius, radius}, {yp, -radius, radius}, {zp, -radius, radius},
		WorkingPrecision -> 100
	];
	by
];

End[]; (* `Private` *)
Protect["longskindepthewjn`*"];
EndPackage[];
