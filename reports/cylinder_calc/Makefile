WS = wolframscript -f
LATEXMK = latexmk -use-make -g -lualatex

all: figures/cylinderBZalongZUniformZ.jpg

.PHONY: all

figures/cylinderBZalongZUniformZ.jpg: calc/cylinderBZalongZUniformZ.csv scripts/cylinderBZalongZUniformZPlot.wls
	$(WS) scripts/cylinderBZalongZUniformZPlot.wls

calc/cylinderBZalongZUniformZ.csv: scripts/cylinderBZalongZUniformZCalc.wls
	$(WS) scripts/cylinderBZalongZUniformZCalc.wls