(* ::Package:: *)

(* :Title: sphereCurrent *)
(* :Context: sphereCurrent` *)
(* :Author: Deepak *)
(* :Date: 2018-10-01 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: *)
(* :Copyright: (c) 2018 Deepak *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["sphereCurrent`", {"longskindepthewjn`"}];
(* Exported symbols added here with SymbolName::usage *)
(* Calculates current for spheres, assuming a constant magnetic field in the Z direction with magnitude 1 *)

outputDirectory = $ScriptCommandLine[[2]];
sphereOutputFilename = $ScriptCommandLine[[3]];
sphereRadius = ToExpression[$ScriptCommandLine[[4]]];
dx = ToExpression[$ScriptCommandLine[[5]]];


Print[TemplateApply[
	StringTemplate["Will create file `4` in directory `1` with current data for sphere of radius `2`, output with a spacing `3`"],
	{outputDirectory, sphereRadius, dx, sphereOutputFilename}
]];
Print[TemplateApply[
	StringTemplate["Integration requires scale factor of `1` for dV, and will generate `2` points on a cubic grid"],
	{dx^3, (2 * ((sphereRadius + 1)/ dx) + 1)^3}
]];

DataDir = outputDirectory;

(* Use longwavelengthewjn script to generate current function *)
region = Ball[{0,0,0}, sphereRadius];
nHat[x_,y_,z_] := Normalize[{x,y,z}];
aSource[x_,y_,z_] := {0, x, 0};
zeroPointPred[x_,y_,z_] := y == -sphereRadius;
j = unscaledJ[aSource, region, zeroPointPred, nHat];
bz = unscaledBz[j, sphereRadius];
Print["Got bz"];
Print[bz[7, 0, 0]];
(* Sample current on a rectangular grid larger than sphere. If we hit any values outside the metal, coerce j to 0*)
(*
jData = Flatten[
	Table[
		Join[{x, y , z}, Quiet[Check[j[x, y , z], {0, 0, 0}, InterpolatingFunction::femdmval], InterpolatingFunction::femdmval]],
		{x, -sphereRadius - 1, sphereRadius + 1, dx},
		{y, -sphereRadius - 1, sphereRadius + 1, dx},
		{z, -sphereRadius - 1, sphereRadius + 1, dx}
	],
	2
];

Print["Finished calculating table, creating sphereJ.csv..."];

Export[
	FileNameJoin[{DataDir, sphereOutputFilename}],
	jData,
	"Table",
	"FieldSeparators" -> ",",
	"LineSeparators" -> "\n",
	DOSTextFormat -> False
];
*)
EndPackage[];
