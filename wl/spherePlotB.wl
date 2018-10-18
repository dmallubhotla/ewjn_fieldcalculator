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
bZ = unscaledBz[j, sphereRadius];
bY = unscaledBy[j, sphereRadius];
bX = unscaledBx[j, sphereRadius];

Print["Calculating bX along diagonal"];
bXAlongDiagonalValues = Table[
	{x, Quiet[bX[x, 0, x], {NIntegrate::slwcon, NIntegrate::eincr}]},
	{x, 3`100, 10`100, 0.5`100}
];
Print[bXAlongDiagonalValues];
dipolePlotAlongDiagonal = Plot[expectedB[x, 0, x][[1]], {x, 2, 7}, PlotLabel -> "B_x along y=0, x=z"];
bXAlongDiagonalValuesPlot = ListPlot[Legended[bXAlongDiagonalValues, "Numerical Integration"], PlotStyle -> Orange];
Print["Exporting plot along diagonal"];
Export["BxAlongDiagonal.jpg", Show[ dipolePlotAlongDiagonal, bXAlongDiagonalValuesPlot], ImageResolution -> 150, ImageSize-> Automatic];


Print["Calculating bZ along x axis"];

bZAlongXAxisValues = Table[
	{x, Quiet[bZ[x, 0, 0], {NIntegrate::slwcon, NIntegrate::eincr}]},
	{x, 4`100, 11`100, 0.5`100}
];
dipolePlotAlongX = Plot[expectedB[x, 0, 0][[3]], {x, 3, 11}, PlotLabel -> "B_z along x-axis"];
bZAlongXAxisValuesPlot = ListPlot[Legended[bZAlongXAxisValues, "Numerical Integration"], PlotStyle -> Orange];
Print["Exporting x axis plot"];
Export["BzAlongXAxis.jpg", Show[ dipolePlotAlongX, bZAlongXAxisValuesPlot], ImageResolution -> 150, ImageSize-> Automatic];

Print["Calculating bZ along diagonal"];
bZAlongDiagonalValues = Table[
	{x, Quiet[bZ[x, 0, x], {NIntegrate::slwcon, NIntegrate::eincr}]},
	{x, 3`100, 10`100, 0.5`100}
];
Print[bZAlongDiagonalValues];
dipolePlotAlongDiagonal = Plot[expectedB[x, 0, x][[3]], {x, 2, 7}, PlotLabel -> "B_z along y=0, x=z"];
bZAlongDiagonalValuesPlot = ListPlot[Legended[bZAlongDiagonalValues, "Numerical Integration"], PlotStyle -> Orange];
Print["Exporting plot along diagonal"];
Export["BzAlongDiagonal.jpg", Show[ dipolePlotAlongDiagonal, bZAlongDiagonalValuesPlot], ImageResolution -> 150, ImageSize-> Automatic];

Print[expectedB[2, 2, 2]];
Print[Quiet[bX[2, 2, 2], {NIntegrate::slwcon, NIntegrate::eincr}]];
Print[Quiet[bY[2, 2, 2], {NIntegrate::slwcon, NIntegrate::eincr}]];
Print[Quiet[bZ[2, 2, 2], {NIntegrate::slwcon, NIntegrate::eincr}]];
EndPackage[];
