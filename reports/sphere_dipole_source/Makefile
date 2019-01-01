WS = wolframscript -f
LATEXMK = latexmk -use-make -g -lualatex

#BZ_FIGURE_FILENAMES = plot1to5.jpg plot5to8.jpg plot8to10.jpg sphereDipolePercentDifferences.jpg
BZ_FIGURE_FILENAMES = plot1to5.jpg
BZ_ALONG_X_FIGURES = $(addprefix figures/bZAlongX, $(BZ_FIGURE_FILENAMES))
BZ_ALONG_Z_FIGURES = $(addprefix figures/bZAlongZ, $(BZ_FIGURE_FILENAMES))

#BZ_ALONG_X_CALCS_FILENAMES = sphereBZalongXDipoleLowRes.csv sphereBZalongXUniform.csv sphereBZalongXDipoleHighRes.csv
BZ_ALONG_X_CALCS_FILENAMES = sphereBZalongXDipoleLowRes.csv
BZ_ALONG_X_CALCS = $(addprefix calc/, $(BZ_ALONG_X_CALCS_FILENAMES))

all: ../output/sphere_dipole_source.pdf

../output/sphere_dipole_source.pdf: sphere_dipole_source.pdf
	cp sphere_dipole_source.pdf ../output/

clean:
	latexmk -C

.PHONY: all clean

sphere_dipole_source.pdf: sphere_dipole_source.tex $(BZ_ALONG_X_FIGURES) $(BZ_ALONG_Z_FIGURES)
	$(LATEXMK) sphere_dipole_source.tex

# sphere dipole source rules
$(BZ_ALONG_X_FIGURES): scripts/sphereBZAlongXPlots.wls $(BZ_ALONG_X_CALCS)
	$(WS) scripts/sphereBZAlongXPlots.wls


$(BZ_ALONG_X_CALCS): scripts/sphereBZAlongXCalcs.wls
	$(WS) scripts/sphereBZAlongXCalcs.wls

$(BZ_ALONG_Z_FIGURES): scripts/sphereBZAlongZPlots.wls calc/sphereBZDipoleHighRes.csv calc/sphereBZalongZUniform.csv
	$(WS) scripts/sphereBZAlongZPlots.wls

calc/sphereBZalongZUniform.csv: scripts/sphereBZalongZUniformCalc.wls
	$(WS) scripts/sphereBZalongZUniformCalc.wls

calc/sphereBZalongZDipoleHighRes.csv: scripts/sphereBZAlongZDipoleHighResCalc.wls
	$(WS) scripts/sphereBZalongZDipoleHighResCalc.wls