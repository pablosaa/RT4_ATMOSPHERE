# Fortran compiler:
CF = gfortran

# Getting server name:
HOST = $(shell hostname)

# netCDF library paths:
ifeq "$(HOST)" "cyclone.hpc.uib.no"
LIBPATH = $(EBROOTNETCDFMINFORTRAN)/lib
INCPATH = $(EBROOTNETCDFMINFORTRAN)/include
else
ifeq ("$(HOST)", $(filter "$(HOST)","gfi3101963.klientdrift.uib.no" "mysu.yolle"))
LIBPATH = /usr/lib
INCPATH = /usr/include
endif
endif

# Base directory where repository is located:
BASE:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

export BASE
export LIBPATH
export INCPATH

OBJ = $(BASE)/obj
BIN = $(BASE)/bin
SRC = $(BASE)/src

$(BIN)/RS_RT3_run: $(OBJ)/%.o
	$(CF) -std=legacy -fopenmp -o $@ $(OBJ)/RS_rt3.o $(OBJ)/rt3.o $(OBJ)/radtran3.o \
        $(OBJ)/radutil3.o  $(OBJ)/radscat3.o $(OBJ)/radintg3.o $(OBJ)/radmat.o \
	$(OBJ)/scat_utilities.o $(OBJ)/dsd_utilities.o $(OBJ)/Fresnel_surf.o \
	$(OBJ)/mitt_time.o $(OBJ)/RT_ncdf_IO.o -I$(INCPATH) -L$(LIBPATH) -lnetcdff

$(OBJ)/%.o: $(SRC)/*.f $(SRC)/*.F90 $(SRC)/*.c 
	+$(MAKE) -C $(SRC)/

	rm -f $(BASE)/output/tmp/*
	@echo '**** IMPORTANT! ****'
	@echo '** If in HPC? do not forget to run> module load netCDF-Fortran'
	@echo '** before start the simulations.'

clean:
	rm -f $(OBJ)/*.o
	rm -f $(BASE)/output/tmp/*.*

#RS_rt3.o rt3.o radtran3.o radutil3.o radscat3.o radintg3.o radmat.o scat_utilities.o dsd_utilities.o Fresnel_surf.o RT_ncdf_IO.o

