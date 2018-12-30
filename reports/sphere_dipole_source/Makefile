WS = wolframscript -f
LATEXMK = latexmk -use-make -g -lualatex

SPHERE_DIPOLE_SOURCE_FIGURE_FILENAMES = plot1to5.jpg plot5to8.jpg plot8to10.jpg sphereDipolePercentDifferences.jpg
SPHERE_DIPOLE_SOURCE_FIGURES = $(addprefix figures/, $(SPHERE_DIPOLE_SOURCE_FIGURE_FILENAMES))

all: ../output/sphere_dipole_source.pdf

../output/sphere_dipole_source.pdf: sphere_dipole_source.pdf
	cp sphere_dipole_source.pdf ../output/

clean:
	latexmk -C

.PHONY: all clean

sphere_dipole_source.pdf: sphere_dipole_source.tex $(SPHERE_DIPOLE_SOURCE_FIGURES)
	$(LATEXMK) sphere_dipole_source.tex

# sphere dipole source rules
$(SPHERE_DIPOLE_SOURCE_FIGURES): scripts/sphereBPlots.wls
	$(WS) scripts/sphereBPlots.wls
