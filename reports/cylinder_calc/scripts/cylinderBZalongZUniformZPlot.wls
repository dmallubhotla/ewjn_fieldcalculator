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

Print[TemplateApply["Will output figures to directory ``", figuresDirectory]];

(* Generated using: *)
points = Import["calc/cylinderBZalongZUniformZ.csv"];
plot = ListPlot[points, ImageSize -> Large];

Export[figure["cylinderBZalongZUniformZ.jpg"], plot, ImageResolution -> 1200];