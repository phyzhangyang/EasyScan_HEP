all: susyhit feynhiggs higgsbounds gm2calc

##########################################################
# SUSYHIT
##########################################################
SUSYHIT_DIR := ./susyhit
SUSYHIT_TAR := susyhit.tar.gz
SUSYHIT_URL := https://www.itp.kit.edu/~maggie/SUSY-HIT/$(SUSYHIT_TAR)

susyhit:
	@echo ">>> Building SUSYHIT"
	@if [ ! -f $(SUSYHIT_TAR) ]; then \
		echo "Downloading SUSYHIT..."; \
		wget $(SUSYHIT_URL); \
	else \
		echo "$(SUSYHIT_TAR) already exists. Skipping download."; \
	fi
	@if [ ! -d $(SUSYHIT_DIR) ]; then \
		mkdir -p $(SUSYHIT_DIR); \
		tar -xf $(SUSYHIT_TAR) -C $(SUSYHIT_DIR); \
		rm -f $(SUSYHIT_TAR); \
	else \
		echo "$(SUSYHIT_DIR) already exists. Skipping extraction."; \
	fi
	@if [ ! -f $(SUSYHIT_DIR)/run ]; then \
		cd $(SUSYHIT_DIR) && make; \
	fi
	@cp -r for_MSSM_mW.susyhit $(SUSYHIT_DIR)/suspect2_lha.in


##########################################################
# FeynHiggs
##########################################################
FH_DIR := ./FeynHiggs
FH_VERSION := 2.19.0
FH_TAR := FeynHiggs-$(FH_VERSION).tar.gz
FH_URL := https://wwwth.mpp.mpg.de/members/heinemey/feynhiggs/newversion/$(FH_TAR)
FH_LIB := $(FH_DIR)/$(shell uname -m)-$(shell uname -s)/lib64/libFH.a

feynhiggs:
	@echo ">>> Building FeynHiggs"
	@if [ ! -f $(FH_TAR) ]; then \
		echo "Downloading FeynHiggs..."; \
		wget --no-check-certificate $(FH_URL); \
	else \
		echo "$(FH_TAR) already exists. Skipping download."; \
	fi
	@if [ ! -d $(FH_DIR) ]; then \
		tar -xf $(FH_TAR); \
		mv FeynHiggs-$(FH_VERSION) $(FH_DIR); \
		rm -f $(FH_TAR); \
	else \
		echo "$(FH_DIR) already exists. Skipping extraction."; \
	fi
	@if [ ! -f $(FH_LIB) ]; then \
		cd $(FH_DIR) && ./configure && make && make install; \
	fi


##########################################################
# HiggsBounds
##########################################################
HB_DIR := ./HiggsBounds
HB_VERSION := 5.10.2
HB_TAR := higgsbounds-$(HB_VERSION).tar.gz
HB_URL := https://gitlab.com/higgsbounds/higgsbounds/-/archive/$(HB_VERSION)/$(HB_TAR)
HB_BIN := $(HB_DIR)/build/example_programs/HBSLHAinputblocksfromFH

higgsbounds:
	@echo ">>> Building HiggsBounds"
	@if [ ! -f $(HB_TAR) ]; then \
		echo "Downloading HiggsBounds..."; \
		wget $(HB_URL); \
	else \
		echo "$(HB_TAR) already exists. Skipping download."; \
	fi
	@if [ ! -d $(HB_DIR) ]; then \
		tar -xf $(HB_TAR); \
		mv higgsbounds-$(HB_VERSION) $(HB_DIR); \
		cp hbounds.txt $(HB_DIR)/CMakeLists.txt; \
		rm -f $(HB_TAR); \
	else \
		echo "$(HB_DIR) already exists. Skipping extraction."; \
	fi
	@if [ ! -f $(HB_BIN) ]; then \
		cd $(HB_DIR) && mkdir -p build && cd build && cmake -DFeynHiggs_ROOT=$$(pwd)/../../FeynHiggs .. && make; \
	fi


##########################################################
# GM2Calc
##########################################################
GM2_DIR := ./GM2Calc
GM2_VERSION := 2.3.1
GM2_TAR := v$(GM2_VERSION).tar.gz
GM2_URL := https://github.com/GM2Calc/GM2Calc/archive/$(GM2_TAR)
GM2_BIN := $(GM2_DIR)/bin/gm2calc.x

gm2calc:
	@echo ">>> Building GM2Calc"
	@if [ ! -f $(GM2_TAR) ]; then \
		echo "Downloading GM2Calc..."; \
		wget $(GM2_URL); \
	else \
		echo "$(GM2_TAR) already exists. Skipping download."; \
	fi
	@if [ ! -d $(GM2_DIR) ]; then \
		tar -xf $(GM2_TAR); \
		mv GM2Calc-$(GM2_VERSION) $(GM2_DIR); \
		rm -f $(GM2_TAR); \
	else \
		echo "$(GM2_DIR) already exists. Skipping extraction."; \
	fi
	@if [ ! -f $(GM2_BIN) ]; then \
		cd $(GM2_DIR) && cmake . && make; \
	fi


##########################################################
# Clean
##########################################################
.PHONY: clean
clean:
	rm -rf $(SUSYHIT_DIR) $(FH_DIR) $(HB_DIR) $(GM2_DIR)
