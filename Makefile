WS = wolframscript -f

all: install

# Convenience targets
.PHONY: install repl all installWL



install: installWL

installWL:
	@$(eval MATHEMATICA_INSTALL_LOCATION=$(shell wolframscript -c 'FileNameJoin[{StringReplace[$$UserBaseDirectory, "\\" -> "/"], "Applications", "longskindepthewjn"}, OperatingSystem -> "Unix"]'))
	@echo $(MATHEMATICA_INSTALL_LOCATION)
	mkdir -p $(MATHEMATICA_INSTALL_LOCATION)
	cp src/wl/longskindepthewjn.wl $(MATHEMATICA_INSTALL_LOCATION)