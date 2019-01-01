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

figure[filename_] := FileNameJoin[{figuresDirectory, filename}];

bZAlongZDipolePoints = Import["calc/sphereBZalongZDipoleHighRes.csv"];
bZAlongXDipolePoints = Import["calc/sphereBZalongXDipoleHighRes.csv"];

zxRatio = MapThread[
	{#1[[1]], #1[[2]] / #2[[2]] } &,
	{bZAlongZDipolePoints, bZAlongXDipolePoints}
];

zxRatioPlot = ListPlot[zxRatio, PlotRange -> {{2, 10}, {0, 6}}, ImageSize -> Large];
Export[figure["sphereDipoleSourceZXRatio.jpg"], zxRatioPlot, ImageResolution -> 1200];