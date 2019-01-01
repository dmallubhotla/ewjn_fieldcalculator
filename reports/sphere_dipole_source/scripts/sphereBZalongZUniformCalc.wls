Needs["longskindepthewjn`"];

sphereUniformbZAlongZPoints = ParallelTable[{z, sphereUniformField[1, 250, 2/(z^3)]["bZ"][0, 0, z]}, {z, 2, 10, .5}];
Export[calc["sphereBZalongZUniform.csv"], sphereUniformbZAlongZPoints];
Print["Exported uniform calc"];