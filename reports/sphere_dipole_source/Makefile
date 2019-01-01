WS = wolframscript -f
LATEXMK = latexmk -use-make -g -lualatex

BZ_ALONG_X_FIGURE_FILENAMES = plot1to5.jpg plot5to8.jpg plot8to10.jpg sphereDipolePercentDifferences.jpg
BZ_ALONG_X_FIGURES = $(addprefix figures/bZAlongX, $(BZ_ALONG_X_FIGURES_FILENAMES))

BZ_ALONG_X_CALCS_FILENAMES = sphereBZalongXDipoleLowRes.csv sphereBZalongXUniform.csv sphereBZalongXDipoleHighRes.csv
BZ_ALONG_X_CALCS = $(addprefix calc/, $(BZ_ALONG_X_CALCS_FILENAMES))

all: ../output/sphere_dipole_source.pdf

../output/sphere_dipole_source.pdf: sphere_dipole_source.pdf
	cp sphere_dipole_source.pdf ../output/

clean:
	latexmk -C

.PHONY: all clean

sphere_dipole_source.pdf: sphere_dipole_source.tex $(BZ_ALONG_X_FIGURES)
	$(LATEXMK) sphere_dipole_source.tex

# sphere dipole source rules
$(SPHERE_DIPOLE_SOURCE_FIGURES): scripts/sphereBZAlongXPlots.wls $(BZ_ALONG_X_CALCS)
	$(WS) scripts/sphereBZAlongXPlots.wls

$(BZ_ALONG_X_CALCS): scripts/sphereBZAlongXCalcs.wls
	$(WS) scripts/sphereBZAlongXCalcs.wls