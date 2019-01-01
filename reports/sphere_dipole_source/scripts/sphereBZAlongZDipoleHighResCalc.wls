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

Print[TemplateApply["Will output calculated values to directory ``", calculatedDirectory]];


sphereDipoleHighResbZAlongZPoints = ParallelTable[{z, sphereDipoleFieldHigherMeshResolution[1, 250, {0, 0, z}, {0, 0 ,1}]["bZ"][0, 0, z]}, {z, 2, 10, .5}];
Export[calc["sphereBZalongZDipoleHighRes.csv"], sphereDipoleHighResbZAlongZPoints];
Print["Exported high res dipole calc"];
