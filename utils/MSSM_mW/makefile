SHELL := /bin/bash

MAKE ?= make
CMAKE ?= cmake
BREW ?= $(shell command -v brew 2>/dev/null)

BOOST_PREFIX ?= $(shell if [ -n "$(BREW)" ]; then "$(BREW)" --prefix boost 2>/dev/null; fi)
EIGEN3_PREFIX ?= $(shell if [ -n "$(BREW)" ]; then "$(BREW)" --prefix eigen@3 2>/dev/null; fi)
EIGEN_PREFIX ?= $(shell if [ -n "$(BREW)" ]; then "$(BREW)" --prefix eigen 2>/dev/null; fi)

GM2_EIGEN3_DIR ?= $(firstword $(wildcard $(EIGEN3_PREFIX)/share/eigen3/cmake $(EIGEN_PREFIX)/share/eigen3/cmake))
GM2_CMAKE_ARGS := -DENABLE_MATHEMATICA=OFF -DENABLE_PYTHON=OFF -DENABLE_TESTS=OFF -DENABLE_EXAMPLES=OFF
ifneq ($(BOOST_PREFIX),)
GM2_CMAKE_ARGS += -DBOOST_ROOT=$(BOOST_PREFIX)
endif
ifneq ($(GM2_EIGEN3_DIR),)
GM2_CMAKE_ARGS += -DEigen3_DIR=$(GM2_EIGEN3_DIR)
endif
GM2_CMAKE_EXTRA_ARGS ?=

.PHONY: all check-env susyhit feynhiggs higgsbounds gm2calc clean

all: check-env susyhit feynhiggs higgsbounds gm2calc

define download_file
if command -v curl >/dev/null 2>&1; then \
	curl -L -o "$(1)" "$(2)"; \
elif command -v wget >/dev/null 2>&1; then \
	wget -O "$(1)" "$(2)"; \
else \
	echo "ERROR: neither curl nor wget was found. Please install one of them or place $(1) in this directory."; \
	exit 1; \
fi
endef

check-env:
	@echo ">>> Checking build environment"
	@command -v "$(MAKE)" >/dev/null 2>&1 || { echo "ERROR: make was not found."; exit 1; }
	@command -v "$(CMAKE)" >/dev/null 2>&1 || { echo "ERROR: cmake was not found."; exit 1; }
	@command -v gfortran >/dev/null 2>&1 || { echo "ERROR: gfortran was not found."; exit 1; }
	@if ! command -v curl >/dev/null 2>&1 && ! command -v wget >/dev/null 2>&1; then \
		echo "ERROR: neither curl nor wget was found. Please install one of them."; \
		exit 1; \
	fi
	@if command -v file >/dev/null 2>&1; then \
		host_arch="$$(uname -m)"; \
		gfortran_path="$$(command -v gfortran)"; \
		gfortran_file="$$(file "$$gfortran_path")"; \
		case "$$host_arch:$$gfortran_file" in \
			arm64:*arm64*|arm64:*aarch64*) ;; \
			aarch64:*aarch64*|aarch64:*arm64*) ;; \
			x86_64:*x86_64*|x86_64:*x86-64*) ;; \
			*) \
				echo "ERROR: gfortran architecture does not match the current shell architecture."; \
				echo "       shell architecture: $$host_arch"; \
				echo "       gfortran: $$gfortran_file"; \
				echo "       On Apple Silicon, install arm64 gcc/gfortran with /opt/homebrew,"; \
				echo "       or run the whole build in Rosetta with: arch -x86_64 zsh"; \
				exit 1; \
				;; \
		esac; \
	fi

##########################################################
# SUSYHIT
##########################################################
SUSYHIT_DIR := ./susyhit
SUSYHIT_TAR := susyhit.tar.gz
SUSYHIT_URL := https://www.itp.kit.edu/~maggie/SUSY-HIT/$(SUSYHIT_TAR)

susyhit:
	@echo ">>> Building SUSYHIT"
	@if [ -d $(SUSYHIT_DIR) ]; then \
		echo "$(SUSYHIT_DIR) already exists. Skipping download."; \
	elif [ ! -f $(SUSYHIT_TAR) ]; then \
		echo "Downloading SUSYHIT..."; \
		$(call download_file,$(SUSYHIT_TAR),$(SUSYHIT_URL)); \
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
		$(MAKE) -C $(SUSYHIT_DIR); \
	else \
		echo "$(SUSYHIT_DIR)/run already exists. Skipping build."; \
	fi
	@cp -r for_MSSM_mW.susyhit $(SUSYHIT_DIR)/suspect2_lha.in


##########################################################
# FeynHiggs
##########################################################
FH_DIR := ./FeynHiggs
FH_VERSION := 2.19.0
FH_TAR := FeynHiggs-$(FH_VERSION).tar.gz
FH_URL := https://wwwth.mpp.mpg.de/members/heinemey/feynhiggs/newversion/$(FH_TAR)
FH_ARCH_DIR := $(FH_DIR)/$(shell uname -m)-$(shell uname -s)
FH_LIB := $(firstword $(wildcard $(FH_ARCH_DIR)/lib/libFH.a $(FH_ARCH_DIR)/lib64/libFH.a))
FH_LIB_DISPLAY := $(if $(FH_LIB),$(FH_LIB),$(FH_ARCH_DIR)/lib/libFH.a or $(FH_ARCH_DIR)/lib64/libFH.a)

feynhiggs:
	@echo ">>> Building FeynHiggs"
	@if [ -d $(FH_DIR) ]; then \
		echo "$(FH_DIR) already exists. Skipping download."; \
	elif [ ! -f $(FH_TAR) ]; then \
		echo "Downloading FeynHiggs..."; \
		$(call download_file,$(FH_TAR),$(FH_URL)); \
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
	@if [ -z "$(FH_LIB)" ]; then \
		cd $(FH_DIR) && ./configure && $(MAKE) && $(MAKE) install; \
	else \
		echo "$(FH_LIB_DISPLAY) already exists. Skipping build."; \
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
	@if [ -d $(HB_DIR) ]; then \
		echo "$(HB_DIR) already exists. Skipping download."; \
	elif [ ! -f $(HB_TAR) ]; then \
		echo "Downloading HiggsBounds..."; \
		$(call download_file,$(HB_TAR),$(HB_URL)); \
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
		cp hbounds.txt $(HB_DIR)/CMakeLists.txt; \
		cd $(HB_DIR) && mkdir -p build && cd build && $(CMAKE) -DFeynHiggs_ROOT=$$(pwd)/../../FeynHiggs .. && $(MAKE); \
	else \
		echo "$(HB_BIN) already exists. Skipping build."; \
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
	@if [ -d $(GM2_DIR) ]; then \
		echo "$(GM2_DIR) already exists. Skipping download."; \
	elif [ ! -f $(GM2_TAR) ]; then \
		echo "Downloading GM2Calc..."; \
		$(call download_file,$(GM2_TAR),$(GM2_URL)); \
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
		if [ -z "$(EIGEN3_PREFIX)" ] && [ -n "$(EIGEN_PREFIX)" ]; then \
			echo "eigen@3 was not found; allowing GM2Calc to use the installed Eigen package."; \
			perl -0pi -e 's/find_package\(Eigen3 3\.1 REQUIRED NO_MODULE\)/find_package(Eigen3 REQUIRED NO_MODULE)/' $(GM2_DIR)/CMakeLists.txt; \
		fi; \
		cd $(GM2_DIR) && $(CMAKE) $(GM2_CMAKE_ARGS) $(GM2_CMAKE_EXTRA_ARGS) . && $(MAKE); \
	else \
		echo "$(GM2_BIN) already exists. Skipping build."; \
	fi


##########################################################
# Clean
##########################################################
clean:
	rm -rf $(SUSYHIT_DIR) $(FH_DIR) $(HB_DIR) $(GM2_DIR)
