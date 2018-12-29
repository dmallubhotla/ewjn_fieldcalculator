WS = wolframscript -f

# Sphere dipole source reports
SPHERE_DIPOLE_SOURCE_FIGURE_FILENAMES = plot1to5.jpg plot5to8.jpg plot8to10.jpg sphereDipolePercentDifferences.jpg
SPHERE_DIPOLE_SOURCE_FIGURES = $(addprefix reports/sphere_dipole_source/figures/, $(SPHERE_DIPOLE_SOURCE_FIGURE_FILENAMES))

all: $(SPHERE_DIPOLE_SOURCE_FIGURES)
	@echo "Done."

# Convenience targets
.PHONY: install repl all installWL



install: installWL

installWL:
	@$(eval MATHEMATICA_INSTALL_LOCATION=$(shell wolframscript -c 'FileNameJoin[{StringReplace[$$UserBaseDirectory, "\\" -> "/"], "Applications", "longskindepthewjn"}, OperatingSystem -> "Unix"]'))
	@echo $(MATHEMATICA_INSTALL_LOCATION)
	mkdir -p $(MATHEMATICA_INSTALL_LOCATION)
	cp src/wl/longskindepthewjn.wl $(MATHEMATICA_INSTALL_LOCATION)


# sphere dipole source rules
$(SPHERE_DIPOLE_SOURCE_FIGURES): reports/sphere_dipole_source/scripts/sphereBPlots.wls
	$(WS) reports/sphere_dipole_source/scripts/sphereBPlots.wls