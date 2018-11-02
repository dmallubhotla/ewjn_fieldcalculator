(* ::Package:: *)

(* :Title: sphereB *)
(* :Context: sphereB` *)
(* :Author: Deepak *)
(* :Date: 2018-10-01 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: *)
(* :Copyright: (c) 2018 Deepak *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["sphereB`", {"longskindepthewjn`"}];
(* Plots calculated B fields *)

expectedB[x_, y_, z_] :=
	(2* Pi*(3.^5)/(15)) *
	((3 * Normalize[{x, y, z}] * Dot[Normalize[{x, y, z}], {0, 0, 1} ]) - {0, 0, 1}) / (Norm[{x, y, z}]^3);

sphereRadius = 3;
region = Ball[{0,0,0}, sphereRadius];
nHat[x_,y_,z_] := Normalize[{x,y,z}];
aSource[x_,y_,z_] := {0, x, 0};
zeroPointPred[x_,y_,z_] := y == -sphereRadius;
j = unscaledJ[aSource, region, zeroPointPred, nHat];
bZ = unscaledBz[j];
bY = unscaledBy[j];
bX = unscaledBx[j];

Print["Calculating bX along diagonal"];
bXAlongDiagonalValues = Table[
	{x, bX[x, 0, x]},
	{x, 3, 10, 1}
];

bXRelativeErrorsAlongDiagonal = Map[
	{#[[1]], Abs[(expectedB[#[[1]], 0, #[[1]]][[1]] - #[[2]]) / expectedB[#[[1]], 0, #[[1]]][[1]]]} &, bXAlongDiagonalValues
];
dipolePlotAlongDiagonal = Plot[expectedB[x, 0, x][[1]], {x, 2, 7}, PlotLabel -> "B_x along y=0, x=z"];
bXAlongDiagonalValuesPlot = ListPlot[Legended[bXAlongDiagonalValues, "Numerical Integration"], PlotStyle -> Orange];

bXAlongDiagonalValuesErrorPlot = ListPlot[Legended[bXRelativeErrorsAlongDiagonal, "New method"], PlotStyle -> Orange];
Print["Exporting plot along diagonal"];
Export["BxAlongDiagonalErrorsComparison.jpg", Show[ bXAlongDiagonalValuesErrorPlot], ImageResolution -> 150, ImageSize-> Automatic];
Export["BxAlongDiagonal.jpg", Show[dipolePlotAlongDiagonal, bXAlongDiagonalValuesPlot], ImageResolution -> 150, ImageSize-> Automatic];

Print["Calculating bZ along x axis"];

bZAlongXAxisValues = Table[
	{x, bZ[x, 0, 0]},
	{x, 4, 11, 1/2}
];
dipolePlotAlongX = Plot[expectedB[x, 0, 0][[3]], {x, 3, 11}, PlotLabel -> "B_z along x-axis"];
bZAlongXAxisValuesPlot = ListPlot[Legended[bZAlongXAxisValues, "Numerical Integration"], PlotStyle -> Orange];
Print["Exporting x axis plot"];
Export["BzAlongXAxis.jpg", Show[ dipolePlotAlongX, bZAlongXAxisValuesPlot], ImageResolution -> 150, ImageSize-> Automatic];

Print["Calculating bZ along diagonal"];
bZAlongDiagonalValues = Table[
	{x, bZ[x, 0, x]},
	{x, 3, 10, 1/2}
];
dipolePlotAlongDiagonal = Plot[expectedB[x, 0, x][[3]], {x, 2, 7}, PlotLabel -> "B_z along y=0, x=z"];
bZAlongDiagonalValuesPlot = ListPlot[Legended[bZAlongDiagonalValues, "Numerical Integration"], PlotStyle -> Orange];
Print["Exporting plot along diagonal"];
Export["BzAlongDiagonal.jpg", Show[ dipolePlotAlongDiagonal, bZAlongDiagonalValuesPlot], ImageResolution -> 150, ImageSize-> Automatic];

EndPackage[];
