Needs["longskindepthewjn`"];

(* Populates the figures directory as ../figures *)
calculatedDirectory = FileNameJoin[{
	ParentDirectory[
		DirectoryName[
			FileNameJoin[{
				Directory[],
				$ScriptCommandLine[[1]]
			}]
		]
	], "calc"
}];

calc[filename_] := FileNameJoin[{calculatedDirectory, filename}];

cylinderUniformbZAlongZPoints = ParallelTable[{z, cylinderUniformField[3, 6, 250, 2/(z^3)]["bZ"][0, 0, z]}, {z, 2, 10, .5}];
Export[calc["cylinderBZalongZUniformZ.csv"], cylinderUniformbZAlongZPoints];
Print["Exported uniform calc"];