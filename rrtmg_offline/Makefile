FC = gfortran

ECRAD_DIR = $(abspath ../ecrad_work/ecrad-1.5.0_source)
ECRAD_MOD = $(ECRAD_DIR)/mod
ECRAD_LIBS = $(ECRAD_DIR)/lib/libradiation.a $(ECRAD_DIR)/lib/libifsrrtm.a \
             $(ECRAD_DIR)/lib/libutilities.a $(ECRAD_DIR)/lib/libifsaux.a \
             $(ECRAD_DIR)/lib/libdrhook.a

# Em sistemas Cray/HPE (compilador chamado via wrapper "ftn"), os
# includes/libs do NetCDF costumam ser injetados automaticamente pelo
# wrapper a partir dos módulos carregados (module load netcdf...) --
# nf-config normalmente não é necessário nesses sistemas e, se estiver
# mal configurado, pode gerar flags quebradas (ex.: "-I" sem caminho,
# que quebra a linha de compilação inteira). Por isso NÃO chamamos
# nf-config por padrão. Em sistemas onde isso for necessário (ex.: um
# desktop/servidor comum com gfortran "puro"), defina explicitamente:
#   make NETCDF_FFLAGS="$(shell nf-config --fflags)" NETCDF_LIBS="$(shell nf-config --flibs)"
# ou edite estas duas linhas diretamente.
NETCDF_DIR = $(shell nf-config --prefix)
NETCDF_INC = $(shell nf-config --includedir)
NETCDF_FFLAGS =$(shell nf-config --fflags)
NETCDF_LIBS   = -L$(NETCDF_DIR)/lib -lnetcdf -lnetcdff
#FFLAGS_DEBUG =  -fconvert=big-endian -fbounds-check -Waliasing -Wall -fbacktrace -ffpe-trap=invalid,overflow,zero  -finit-real=nan -finit-integer=nan  -ffree-line-length-none  -O0 -Warray-bounds  -ffast-math -funroll-loops -ftree-vectorizer-verbose=2
#FFLAGS_DEBUG =  -fcheck=all -fbounds-check -Waliasing -Wall -fbacktrace  -finit-real=nan -finit-integer=nan  -ffree-line-length-none  -O0 -Warray-bounds  -ffast-math -funroll-loops -ftree-vectorizer-verbose=2

FFLAGS =  -cpp  -finit-real=nan -finit-integer=nan -ffree-form -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -fconvert=big-endian -O2 \
         $(FFLAGS_DEBUG) $(NETCDF_FFLAGS) -I$(NETCDF_INC) -I$(ECRAD_MOD)
LDFLAGS = $(NETCDF_LIBS) -fopenmp

SRC = src
BUILD = build

# Compile order matters (module dependencies). RRTMG's own files first
# (no ecRad dependency), then the ecRad bridge (needs $(ECRAD_MOD) on
# the include path, and $(ECRAD_LIBS) at link time), then the shared
# namelist reader, then the two driver programs.
OBJS = mpas_stubs_offline.o \
       module_ra_rrtmg_vinterp.o \
       module_ra_rrtmg_lw.o \
       module_ra_rrtmg_sw.o \
       mpas_atmphys_rrtmg_lwinit.o \
       mpas_atmphys_rrtmg_swinit.o \
       radconst_offline.o \
       era5_bridge.o \
       era5_ecrad_bridge.o \
       output_netcdf.o \
       namelist_offline.o

.PHONY: all clean test ecrad-lib

all: ecrad-lib $(BUILD)/rad_offline_driver $(BUILD)/rrtmg_offline_driver

# Builds the ecRad library by invoking ecRad's OWN Makefile (this
# project deliberately does not duplicate ecRad's build system, which
# has its own PROFILE-based compiler configuration in
# Makefile_include.<profile> -- see $(ECRAD_DIR)/Makefile). Skips the
# invocation if libradiation.a is already there, so a second `make`
# here doesn't needlessly rebuild ecRad.
#
# ECRAD_PROFILE selects which Makefile_include.<profile> ecRad itself
# uses (gfortran by default); override on the command line, e.g.
# `make ECRAD_PROFILE=ifs` if you build with a different compiler.
ECRAD_PROFILE ?= gfortran

ecrad-lib:
	@if [ ! -f $(ECRAD_DIR)/lib/libradiation.a ]; then \
	  if [ ! -d $(ECRAD_DIR) ]; then \
	    echo "ERROR: $(ECRAD_DIR) not found."; \
	    echo "Extract the ecRad 1.5.0 source there first (see DADOS_ECRAD.md), or set ECRAD_DIR."; \
	    exit 1; \
	  fi; \
	  echo "Applying ecrad_patches/ (compile-bug and ozone-related fixes -- see DADOS_ECRAD.md)..."; \
	  cp ecrad_patches/radiation_interface.F90 $(ECRAD_DIR)/radiation/radiation_interface.F90; \
	  cp ecrad_patches/module_ra_ecrad_sw.F90   $(ECRAD_DIR)/radiation/module_ra_ecrad_sw.F90; \
	  cp ecrad_patches/module_ra_ecrad_lw.F90   $(ECRAD_DIR)/radiation/module_ra_ecrad_lw.F90; \
	  echo "libradiation.a not found -- building ecRad (PROFILE=$(ECRAD_PROFILE)) in $(ECRAD_DIR)..."; \
	  $(MAKE) -C $(ECRAD_DIR) PROFILE=$(ECRAD_PROFILE); \
	fi

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/%.o: $(SRC)/%.F90 | $(BUILD)
	cd $(BUILD) && $(FC) $(FFLAGS) -c ../$<

$(BUILD)/%.o: $(SRC)/%.F | $(BUILD)
	cd $(BUILD) && $(FC) $(FFLAGS) -c ../$<

$(BUILD)/rad_offline_driver: $(addprefix $(BUILD)/,$(OBJS)) $(BUILD)/rad_offline_driver.o
	cd $(BUILD) && $(FC) -o rad_offline_driver $(OBJS) rad_offline_driver.o \
	  $(ECRAD_LIBS) $(LDFLAGS)

# Kept for anyone who only wants RRTMG and prefers the plain
# command-line interface over a namelist file; does not link against
# ecRad at all.
$(BUILD)/rrtmg_offline_driver: $(addprefix $(BUILD)/,mpas_stubs_offline.o module_ra_rrtmg_vinterp.o \
	  module_ra_rrtmg_lw.o module_ra_rrtmg_sw.o mpas_atmphys_rrtmg_lwinit.o \
	  mpas_atmphys_rrtmg_swinit.o radconst_offline.o era5_bridge.o) \
	  $(BUILD)/rrtmg_offline_driver.o
	cd $(BUILD) && $(FC) -o rrtmg_offline_driver mpas_stubs_offline.o \
	  module_ra_rrtmg_vinterp.o module_ra_rrtmg_lw.o module_ra_rrtmg_sw.o \
	  mpas_atmphys_rrtmg_lwinit.o mpas_atmphys_rrtmg_swinit.o \
	  radconst_offline.o era5_bridge.o rrtmg_offline_driver.o $(LDFLAGS)

test: $(addprefix $(BUILD)/,mpas_stubs_offline.o era5_bridge.o)
	cd $(BUILD) && $(FC) $(FFLAGS) -c ../$(SRC)/test_era5_bridge.F90
	cd $(BUILD) && $(FC) -o test_era5_bridge mpas_stubs_offline.o era5_bridge.o test_era5_bridge.o $(LDFLAGS)
	cd $(BUILD) && ./test_era5_bridge

clean:
	 @echo "$(MAKE)"
	 rm -rf $(BUILD)     
	 (cd $(ECRAD_DIR);$(MAKE) -f Makefile clean)
