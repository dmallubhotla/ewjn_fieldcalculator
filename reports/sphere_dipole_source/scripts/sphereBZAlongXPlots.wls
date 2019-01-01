(* Populates the figures directory as ../figures *)
figuresDirectory = FileNameJoin[{
	ParentDirectory[
		DirectoryName[
			FileNameJoin[{
				Directory[],
				$ScriptCommandLine[[1]]
			}]
		]
	], "figures"
}];

figure[filename_] := FileNameJoin[{figuresDirectory, "bZAlongX" <> filename}];

Print[TemplateApply["Will output figures to directory ``", figuresDirectory]];

(* Generated as follows: *)
(* Table[{x, sphereDipoleField[1, 250, {x, 0, 0}, {0, 0 ,1}]["bZ"][x, 0, 0]}, {x, 2, 10, .5}] *)
sphereDipoleFieldPoints = Import["calc/sphereBZAlongXDipoleLowRes.csv"];

(* Generated using: *)
(*sphereUniformFieldPoints = Table[{x, sphereUniformField[1, 250, 1/(x^3)]["bZ"][x, 0, 0]}, {x, 2, 10, .5}]*)
sphereUniformFieldPoints = Import["calc/sphereBZAlongXUniform.csv"];

(* ParallelTable[{x, sphereDipoleFieldHigherMeshResolution[1, 250, {x, 0, 0}, {0, 0, 1}]["bZ"][x, 0, 0]}, {x, 2.5, 10, .5}]*)
sphereDipoleHigherMeshResolutionPoints = Import["calc/sphereBZAlongXDipoleHighRes.csv"];

expectedBz[x_] := 1/(15 x^6 * 250^2);

Print["Generating plots..."];

plot1to5 = Show[Plot[expectedBz[x], {x, 1.5, 5} (*, PlotLabels \[Rule]
  "Exact uniform field sphere" *)],
	ListPlot[<|
		"Uniform source" -> sphereUniformFieldPoints,
		"Dipole source, low mesh resolution" -> sphereDipoleFieldPoints,
		"Dipole source, high mesh resolution" -> sphereDipoleHigherMeshResolutionPoints |>], ImageSize -> Large];

Export[figure["plot1to5.jpg"], plot1to5, ImageResolution -> 1200];

plot5to8 = Show[Plot[expectedBz[x], {x, 5, 8} (*, PlotLabels \[Rule]
  "Exact uniform field sphere" *)],
	ListPlot[<|
		"Uniform source" -> sphereUniformFieldPoints,
		"Dipole source, low mesh resolution" -> sphereDipoleFieldPoints,
		"Dipole source, high mesh resolution" -> sphereDipoleHigherMeshResolutionPoints |>], ImageSize -> Large];

Export[figure["plot5to8.jpg"], plot5to8, ImageResolution -> 1200];

plot8to10 = Show[Plot[expectedBz[x], {x, 8, 10} (*, PlotLabels \[Rule]
  "Exact uniform field sphere" *)],
	ListPlot[<|
		"Uniform source" -> sphereUniformFieldPoints,
		"Dipole source, low mesh resolution" -> sphereDipoleFieldPoints,
		"Dipole source, high mesh resolution" -> sphereDipoleHigherMeshResolutionPoints |>], ImageSize -> Large];


Export[figure["plot8to10.jpg"], plot8to10, ImageResolution -> 1200];

Print["Now generating error plots..."];

lowResolutionErrors = MapThread[{#1[[
		1]], (#2[[2]] - #1[[2]]) / (#2[[2]]) } &, {sphereDipoleFieldPoints,
	sphereUniformFieldPoints}];

highResolutionErrors = MapThread[{#1[[
		1]], (#2[[2]] - #1[[2]]) / (#2[[2]]) } &, {sphereDipoleHigherMeshResolutionPoints,
	sphereUniformFieldPoints}];

errorPlot = ListPlot[<| "low mesh resolution percent differences" -> lowResolutionErrors, "high mesh resolution percent differences" -> highResolutionErrors |>, ImageSize -> Large];
Export[figure["sphereDipolePercentDifferences.jpg"], errorPlot, ImageResolution -> 1200];