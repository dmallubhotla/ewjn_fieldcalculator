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

(* Generated as follows: *)
(* Table[{x, sphereDipoleField[1, 250, {x, 0, 0}, {0, 0 ,1}]["bZ"][x, 0, 0]}, {x, 2, 10, .5}] *)
sphereDipoleFieldPoints = {{2.`, 2.6142124585971176`*^-8}, {2.5`,
	5.6198233371519495`*^-9}, {3.`, 1.6850866278381694`*^-9}, {3.5`,
	6.201294361210435`*^-10}, {4.`, 2.625078148376947`*^-10}, {4.5`,
	1.2303620653511336`*^-10}, {5.`, 6.229366903553939`*^-11}, {5.5`,
	3.349431772616367`*^-11}, {6.`, 1.8889472855369908`*^-11}, {6.5`,
	1.106822216445101`*^-11}, {7.`, 6.687668797382179`*^-12}, {7.5`,
	4.140880647713249`*^-12}, {8.`, 2.6131198399413383`*^-12}, {8.5`,
	1.6721669276250677`*^-12}, {9.`, 1.0796309812192286`*^-12}, {9.5`,
	6.995145386687253`*^-13}, {10.`, 4.5191194095435516`*^-13}};

(* Generated using: *)
(*sphereUniformFieldPoints = Table[{x, sphereUniformField[1, 250, 1/(x^3)]["bZ"][x, 0, 0]}, {x, 2, 10, .5}]*)
sphereUniformFieldPoints = {{2.`, 1.649680571469782`*^-8}, {2.5`,
	4.3245129277712425`*^-9}, {3.`, 1.4482687339962462`*^-9}, {3.5`,
	5.743389134522651`*^-10}, {4.`, 2.5776057472289147`*^-10}, {4.5`,
	1.2714559655473036`*^-10}, {5.`, 6.75703801414958`*^-11}, {5.5`,
	3.814171704052086`*^-11}, {6.`, 2.26291826063245`*^-11}, {6.5`,
	1.3998974549172707`*^-11}, {7.`, 8.974042451777574`*^-12}, {7.5`,
	5.932104270870742`*^-12}, {8.`, 4.027508188201787`*^-12}, {8.5`,
	2.7993860729621852`*^-12}, {9.`, 1.9866496888852804`*^-12}, {9.5`,
	1.4362654619010345`*^-12}, {10.`, 1.0557870908750322`*^-12}};

(* ParallelTable[{x, sphereDipoleFieldHigherMeshResolution[1, 250, {x, 0, 0}, {0, 0, 1}]["bZ"][x, 0, 0]}, {x, 2.5, 10, .5}]*)
sphereDipoleHigherMeshResolutionPoints = {{2.`, 2.656181406365209`*^-8}, {2.5`, 5.7524366570620894`*^-9}, {3.`,
	1.741065878674298`*^-9}, {3.5`, 6.480964919680831`*^-10}, {4.`,
	2.7811837796711724`*^-10}, {4.5`, 1.3246565222834076`*^-10}, {5.`,
	6.833788567732815`*^-11}, {5.5`, 3.755270968053709`*^-11}, {6.`,
	2.1718255525669784`*^-11}, {6.5`, 1.3101593781215035`*^-11}, {7.`,
	8.187539010125132`*^-12}, {7.5`, 5.271842089803379`*^-12}, {8.`,
	3.4822526923678008`*^-12}, {8.5`, 2.3512340794555217`*^-12}, {9.`,
	1.61798091795782`*^-12}, {9.5`, 1.1318566440898594`*^-12}, {10.`,
	8.031509314684026`*^-13}};

expectedBz[x_] := 1/(15 x^6 * 250^2);

Print["Generating plots..."];

plot1to5 = Show[Plot[expectedBz[x], {x, 1.5, 5} (*, PlotLabels \[Rule]
  "Exact uniform field sphere" *)],
	ListPlot[<|"Dipole source" -> sphereDipoleFieldPoints,
		"Uniform source" -> sphereUniformFieldPoints,
		"Dipole higher resolution source" -> sphereDipoleHigherMeshResolutionPoints |>], ImageSize -> Large];

Export[figure["plot1to5.jpg"], plot1to5];

plot5to8 = Show[Plot[expectedBz[x], {x, 5, 8} (*, PlotLabels \[Rule]
  "Exact uniform field sphere" *)],
	ListPlot[<|"Dipole source" -> sphereDipoleFieldPoints,
		"Uniform source" -> sphereUniformFieldPoints,
		"Dipole higher resolution source" -> sphereDipoleHigherMeshResolutionPoints |>], ImageSize -> Large];

Export[figure["plot5to8.jpg"], plot5to8];

plot8to10 = Show[Plot[expectedBz[x], {x, 8, 10} (*, PlotLabels \[Rule]
  "Exact uniform field sphere" *)],
	ListPlot[<|"Dipole source" -> sphereDipoleFieldPoints,
		"Uniform source" -> sphereUniformFieldPoints,
		"Dipole higher resolution source" -> sphereDipoleHigherMeshResolutionPoints |>], ImageSize -> Large];


Export[figure["plot8to10.jpg"], plot8to10];

Print["Now generating error plots..."];

lowResolutionErrors = MapThread[{#1[[
		1]], (#2[[2]] - #1[[2]]) / (#2[[2]]) } &, {sphereDipoleFieldPoints,
	sphereUniformFieldPoints}];

highResolutionErrors = MapThread[{#1[[
		1]], (#2[[2]] - #1[[2]]) / (#2[[2]]) } &, {sphereDipoleHigherMeshResolutionPoints,
	sphereUniformFieldPoints}];

errorPlot = ListPlot[<| "low mesh resolution percent differences" -> lowResolutionErrors, "high mesh resolution percent differences" -> highResolutionErrors |>, ImageSize -> Large];
Export[figure["sphereDipolePercentDifferences.jpg"], errorPlot];