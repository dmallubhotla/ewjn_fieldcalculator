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


unscaledJ::usage = "unscaledJ[aSource, region, zeroPoint, nHat] returns { jFunc: J[x,y,z], mesh: ElementMesh } for some 3D region region without prefactors. Tricky parts for the arguments are described in comments";

unscaledBz::usage = "unscaledBz[j] returns a function returning the magnetic field without prefactors in the z direction for a current j as an association containing jFunc and mesh";
unscaledBy::usage = "unscaledBy[j] returns a function returning the magnetic field without prefactors in the y direction for a current j as an association containing jFunc and mesh";
unscaledBx::usage = "unscaledBx[j] returns a function returning the magnetic field without prefactors in the x direction for a current j as an association containing jFunc and mesh";

(* For convenience, these functions, when called, return a source potential function. *)

getDipoleSourcePotential::usage = "getDipoleSourcePotential[{mx_, my_, mz_}, {rx_, ry_, rz_}] returns a function a[x, y, z] returning a triple representing A from a point dipole (mx, my, mz) source at (rx, ry, rz)";
getUniformSourcePotential::usage = "getUniformSourcePotential[] returns a function a[x, y, z] returning a function a[x, y, z] representing a potential for a uniform magnetic field in the -z direction";

(* Wrappers *)
sphereUniformField::usage = "sphereUniformFieldCalc[sphereRadius_, skinDepth_, externalFieldStrength_] returns a calculation context for a sphere of a particular radius with a particular skin depth in a uniform field orientated in the -z direction.";
sphereDipoleField::usage = "sphereUniformFieldCalc[sphereRadius_, skinDepth_, dipolePosition_, dipoleMoment_] returns a calculation context for a sphere of a particular radius with a particular skin depth in a dipole field.";
sphereDipoleFieldHigherMeshResolution::usage = "Higher integration order for j calcs version of sphereDipoleField"
cylinderUniformField::usage = "cylinderUniformField[cylinderRadius_, cylinderHeight_, skinDepth_, externalFieldStrength_]";
cylinderDipoleField::usage = "cylinderDipoleField[cylinderRadius_, cylinderHeight_, skinDepth_, dipolePosition_, {mx_, my_, mz_}]";

Begin["`Private`"];

(*
	zeroPointPred is arbitrary, but it makes Mathematica happier to be able to apply a Dirichlet condition at a single point.
 It must be some sort of function that identifies a single point on the surface. Because the mesh is approximate,
 it's best to use an inequality on a single dimension. For example, for a sphere or cylinder you can specify x = y = 0, z <= 0
 to pick the intersection of the z-axis with the lower part of the surface. This is shape dependent.
 *)
(* nHat[x,y,z] specifies the outward facing surface normal at a point (x, y, z) on the surface.

	Do not use Piecewise, as it will interfere with the NDSolveValue. Use HeavisideTheta instead.

	The values of nHat off the surface don't matter; it will only be queried on the 2D surface. However, the mesh is
	not exact on the surface, so avoid having any sharp transitions near the surface, as they can be evaluated unpredictably.
	HeavisideTheta should be defined slightly inside or slightly outside the surface.
	(That means that surface-specific effects or corners need to be very carefully specified!)
*)
solveForF3D[aSource_, region_, zeroPointPred_, nHat_] := NDSolveValue[{
	Laplacian[f[x, y, z], {x, y, z}] == NeumannValue[Dot[aSource[x, y, z], nHat[x, y, z]], True],
	DirichletCondition[f[x, y, z] == 0, zeroPointPred[x, y, z]]
}, f, Element[{x, y, z}, region]];

solveForF3DHIO[aSource_, region_, zeroPointPred_, nHat_, maxCellMeasure_ : 10000] := NDSolveValue[{
	Laplacian[f[x, y, z], {x, y, z}] == NeumannValue[Dot[aSource[x, y, z], nHat[x, y, z]], True],
	DirichletCondition[f[x, y, z] == 0, zeroPointPred[x, y, z]]
}, f, Element[{x, y, z}, region], Method -> {
		"FiniteElement",
		"IntegrationOrder" -> 5,
		"MeshOptions" -> {
			"MaxCellMeasure" -> 1/maxCellMeasure
		}
	}
];

(* Utility to simplify some of the type handling here, to make sure we can move the piecewise inside the list for j*)
(* This essentially allows us to force j to be zero outside the surface, which is important for plotting and integrating.
 This is a bit of a hack, but shouldn't affect performance. This is particularly important if integrating over regions outside the region. *)
piecewise3DFunc[function_, region_, default_, index_] := Function[{x, y, z},
		Simplify`PWToUnitStep[Piecewise[
			{{function[x, y, z][[index]], {x, y, z} \[Element] region}}
			, default
		]]
	];

unscaledJ[aSource_, region_, zeroPoint_, nHat_] := Module[{f, mesh, gradFPlusAs, j},
	f = solveForF3D[aSource, region, zeroPoint, nHat];
	mesh = f["ElementMesh"];
	gradFPlusAs[x_, y_, z_] := Evaluate[D[f[x, y, z], {{x, y, z}}] + aSource[x, y, z]];
	j[x_, y_, z_] = Through[
			Thread[Unevaluated[piecewise3DFunc[gradFPlusAs, region, {0, 0, 0}, {1, 2, 3}]]]
			[x, y, z]
		];
	<| "jFunc" -> j, "mesh" -> mesh |>
];

unscaledJHIO[aSource_, region_, zeroPoint_, nHat_, maxCellMeasure_ : 10000] := Module[{f, mesh, gradFPlusAs, j},
	f = solveForF3DHIO[aSource, region, zeroPoint, nHat, maxCellMeasure];
	mesh = f["ElementMesh"];
	gradFPlusAs[x_, y_, z_] := Evaluate[D[f[x, y, z], {{x, y, z}}] + aSource[x, y, z]];
	j[x_, y_, z_] = Through[
		Thread[Unevaluated[piecewise3DFunc[gradFPlusAs, region, {0, 0, 0}, {1, 2, 3}]]]
		[x, y, z]
	];
	<| "jFunc" -> j, "mesh" -> mesh |>
];

unscaledBz[j_] := Module[{jf, mesh, bz},
	jf = j["jFunc"];
	mesh = j["mesh"];
	bz[x_, y_, z_] := NIntegrate[
		((jf[xp, yp, zp][[1]]*(y - yp)) - (jf[xp, yp, zp][[2]]*(x - xp)) )/(Norm[{x, y, z} - {xp, yp, zp}]^3),
		Element[{xp, yp, zp}, mesh]
	];
	bz
];

unscaledBx[j_] := Module[{jf, mesh, bz},
	jf = j["jFunc"];
	mesh = j["mesh"];
	bx[x_, y_, z_] := NIntegrate[
		((jf[xp, yp, zp][[2]]*(z - zp)) - (jf[xp, yp, zp][[3]]*(y - yp)) )/(Norm[{x, y, z} - {xp, yp, zp}]^3),
		Element[{xp, yp, zp}, mesh]
	];
	bx
];

unscaledBy[j_] := Module[{jf, mesh, bz},
	jf = j["jFunc"];
	mesh = j["mesh"];
	by[x_, y_, z_] := NIntegrate[
		((jf[xp, yp, zp][[3]]*(x - xp)) - (jf[xp, yp, zp][[1]]*(z - zp)) )/(Norm[{x, y, z} - {xp, yp, zp}]^3),
		Element[{xp, yp, zp}, mesh]
	];
	by
];

sphereUniformField[sphereRadius_, skinDepth_, externalFieldStrength_] := Module[{region, nHat, aSource, zpp, j, bX, bY, bZ, uBx, uBy, uBz},
	region = Ball[{0,0,0}, sphereRadius];
	nHat[x_,y_,z_] := Normalize[{x,y,z}];
	aSource[x_,y_,z_] := {0, x, 0};
	zpp[x_, y_, z_] := (x == 0 && y == 0 && z <= 0);
	j = unscaledJ[aSource, region, zpp, nHat];
	uBx = unscaledBx[j];
	uBy = unscaledBy[j];
	uBz = unscaledBz[j];
	bX[x_, y_, z_] := (externalFieldStrength/(2 * Pi * (skinDepth ^2))) * uBx[x, y, z];
	bY[x_, y_, z_] := (externalFieldStrength/(2 * Pi * (skinDepth ^2))) * uBy[x, y, z];
	bZ[x_, y_, z_] := (externalFieldStrength/(2 * Pi * (skinDepth ^2))) * uBz[x, y, z];
	<| "region" -> region, "current" -> j["jFunc"], "integrationMesh" -> j["mesh"], "bX" -> bX, "bY" -> bY, "bZ" -> bZ |>
];

sphereDipoleField[sphereRadius_, skinDepth_, dipolePosition_, {mx_, my_, mz_}] := Module[{region, nHat, aSource, zpp, j, bX, bY, bZ, uBx, uBy, uBz},
	region = Ball[{0,0,0}, sphereRadius];
	nHat[x_,y_,z_] := Normalize[{x,y,z}];
	aSource = getDipoleSourcePotential[{mx, my, mz}, dipolePosition];
	zpp[x_, y_, z_] := (x == 0 && y == 0 && z <= 0);
	j = unscaledJ[aSource, region, zpp, nHat];
	uBx = unscaledBx[j];
	uBy = unscaledBy[j];
	uBz = unscaledBz[j];
	bX[x_, y_, z_] := (1/(2 * Pi * (skinDepth ^2))) * uBx[x, y, z];
	bY[x_, y_, z_] := (1/(2 * Pi * (skinDepth ^2))) * uBy[x, y, z];
	bZ[x_, y_, z_] := (1/(2 * Pi * (skinDepth ^2))) * uBz[x, y, z];
	<| "region" -> region, "current" -> j["jFunc"], "integrationMesh" -> j["mesh"], "bX" -> bX, "bY" -> bY, "bZ" -> bZ |>
];

sphereDipoleFieldHigherMeshResolution[sphereRadius_, skinDepth_, dipolePosition_, {mx_, my_, mz_}, maxCellMeasure_ : 10000] := Module[{region, nHat, aSource, zpp, j, bX, bY, bZ, uBx, uBy, uBz},
	region = Ball[{0,0,0}, sphereRadius];
	nHat[x_,y_,z_] := Normalize[{x,y,z}];
	aSource = getDipoleSourcePotential[{mx, my, mz}, dipolePosition];
	zpp[x_, y_, z_] := (x == 0 && y == 0 && z <= 0);
	j = unscaledJHIO[aSource, region, zpp, nHat, maxCellMeasure];
	uBx = unscaledBx[j];
	uBy = unscaledBy[j];
	uBz = unscaledBz[j];
	bX[x_, y_, z_] := (1/(2 * Pi * (skinDepth ^2))) * uBx[x, y, z];
	bY[x_, y_, z_] := (1/(2 * Pi * (skinDepth ^2))) * uBy[x, y, z];
	bZ[x_, y_, z_] := (1/(2 * Pi * (skinDepth ^2))) * uBz[x, y, z];
	<| "region" -> region, "current" -> j["jFunc"], "integrationMesh" -> j["mesh"], "bX" -> bX, "bY" -> bY, "bZ" -> bZ |>
];

cylinderUniformField[cylinderRadius_, cylinderHeight_, skinDepth_, externalFieldStrength_] := Module[{region, nHat, aSource, zpp, j, bX, bY, bZ, uBx, uBy, uBz},
	region = Cylinder[{{0, 0, -cylinderHeight / 2}, {0, 0, cylinderHeight / 2}}, cylinderRadius];
	nHat[x_, y_, z_] := {
				(x /Sqrt[x^2 + y^2  + (1 - HeavisideTheta[x^2 + y^2 - cylinderRadius^2])])*
						HeavisideTheta[x^2 + y^2 - .99*cylinderRadius^2 ],
				(y /Sqrt[x^2 + y^2 + (1 - HeavisideTheta[x^2 + y^2 - cylinderRadius^2])])*
						HeavisideTheta[x^2 + y^2 - .99* cylinderRadius^2],
				(Sign[z])* HeavisideTheta[.99*cylinderRadius^2 - x^2 - y^2]
			};
	aSource[x_,y_,z_] := {0, x, 0};
	zpp[x_, y_, z_] := (x == 0 && y == 0 && z <= 0);
	j = unscaledJ[aSource, region, zpp, nHat];
	uBx = unscaledBx[j];
	uBy = unscaledBy[j];
	uBz = unscaledBz[j];
	bX[x_, y_, z_] := (externalFieldStrength/(2 * Pi * (skinDepth ^2))) * uBx[x, y, z];
	bY[x_, y_, z_] := (externalFieldStrength/(2 * Pi * (skinDepth ^2))) * uBy[x, y, z];
	bZ[x_, y_, z_] := (externalFieldStrength/(2 * Pi * (skinDepth ^2))) * uBz[x, y, z];
	<| "region" -> region, "current" -> j["jFunc"], "integrationMesh" -> j["mesh"], "bX" -> bX, "bY" -> bY, "bZ" -> bZ |>
];

cylinderDipoleField[cylinderRadius_, cylinderHeight_, skinDepth_, dipolePosition_, {mx_, my_, mz_}] := Module[{region, nHat, aSource, zpp, j, bX, bY, bZ, uBx, uBy, uBz},
	region = Cylinder[{{0, 0, -cylinderHeight / 2}, {0, 0, cylinderHeight / 2}}, cylinderRadius];
	nHat[x_, y_, z_] := {
		(x /Sqrt[x^2 + y^2  + (1 - HeavisideTheta[x^2 + y^2 - cylinderRadius^2])])*
				HeavisideTheta[x^2 + y^2 - .99*cylinderRadius^2 ],
		(y /Sqrt[x^2 + y^2 + (1 - HeavisideTheta[x^2 + y^2 - cylinderRadius^2])])*
				HeavisideTheta[x^2 + y^2 - .99* cylinderRadius^2],
		(Sign[z])* HeavisideTheta[.99*cylinderRadius^2 -  x^2 - y^2]
	};
	aSource = getDipoleSourcePotential[{mx, my, mz}, dipolePosition];
	zpp[x_, y_, z_] := (x == 0 && y == 0 && z <= 0);
	j = unscaledJ[aSource, region, zpp, nHat];
	uBx = unscaledBx[j];
	uBy = unscaledBy[j];
	uBz = unscaledBz[j];
	bX[x_, y_, z_] := (1/(2 * Pi * (skinDepth ^2))) * uBx[x, y, z];
	bY[x_, y_, z_] := (1/(2 * Pi * (skinDepth ^2))) * uBy[x, y, z];
	bZ[x_, y_, z_] := (1/(2 * Pi * (skinDepth ^2))) * uBz[x, y, z];
	<| "region" -> region, "current" -> j["jFunc"], "integrationMesh" -> j["mesh"], "bX" -> bX, "bY" -> bY, "bZ" -> bZ |>
];

getDipoleSourcePotential[{mx_, my_, mz_}, {rx_, ry_, rz_}] :=
		Module[{as},
			as[x_, y_, z_] :=
					Cross[{mx, my,
						mz}, {x, y, z} - {rx, ry,
						rz}]/(Norm[{x, y, z} - {rx, ry, rz}]^3);
			as
		];

getUniformSourcePotential[] := Module[{as},
			as[x_, y_, z_] := {0, x, 0};
			as
		];


End[]; (* `Private` *)
Protect["longskindepthewjn`*"];
EndPackage[];
