all: repl

# Convenience targets
.PHONY: install repl all mathematicaLocation


	
install: installWL
	poetry install

	
repl: $(ENV)
	poetry run python -i

installWL:
	@$(eval MATHEMATICA_INSTALL_LOCATION=$(shell wolframscript -c 'FileNameJoin[{StringReplace[$$UserBaseDirectory, "\\" -> "/"], "Applications", "longskindepthewjn"}, OperatingSystem -> "Unix"]'))
	@echo $(MATHEMATICA_INSTALL_LOCATION)
	mkdir -p $(MATHEMATICA_INSTALL_LOCATION)
	cp wl/longskindepthewjn.wl $(MATHEMATICA_INSTALL_LOCATION)