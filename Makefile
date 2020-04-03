# Fortran compiler:
FC = $(shell pkg-config --variable=fcompiler netcdf-fortran)
CC = $(shell pkg-config --variable=ccompiler netcdf-fortran)

# Getting needed libraries path:
LIBPATH = $(shell pkg-config --libs netcdf-fortran)
INCPATH = $(shell PKG_CONFIG_ALLOW_SYSTEM_CFLAGS=1 pkg-config --cflags netcdf-fortran)


# Getting server name:
HOSTNAME = $(shell hostname)

# netCDF library paths:
ifeq "$(HOSTNAME)" "cyclone.hpc.uib.no"
MSG = 'module load netCDF-Fortran'
#LIBPATH = $(EBROOTNETCDFMINFORTRAN)/lib
#INCPATH = $(EBROOTNETCDFMINFORTRAN)/include
else
ifeq ("$(HOSTNAME)", $(filter "$(HOSTNAME)","gfi3101963.klientdrift.uib.no" "mysu" ))
MSG = 'No module load need it'
#LIBPATH = /usr/lib
#INCPATH = /usr/include
else
MSG = 'module load gnu \n module load openmpi \n module load netcdf-fortran'
endif
endif

# Base directory where repository is located:
BASE:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

export FC
export BASE
export LIBPATH
export INCPATH

OBJ = $(BASE)/obj
BIN = $(BASE)/bin
SRC = $(BASE)/src

FLAGS = -C -mcmodel=medium -std=legacy -Wall -fcheck=all 

$(BIN)/RS_RT3_run: $(OBJ)/%.o
	$(FC) $(FLAGS) -o $@ $(OBJ)/RS_rt3.o $(OBJ)/rt3.o $(OBJ)/radtran3.o \
        $(OBJ)/radutil3.o  $(OBJ)/radscat3.o $(OBJ)/radintg3.o $(OBJ)/radmat.o \
	$(OBJ)/scat_utilities.o $(OBJ)/dsd_utilities.o $(OBJ)/Fresnel_surf.o \
	$(OBJ)/mitt_time.o $(OBJ)/RT_ncdf_IO.o $(INCPATH) $(LIBPATH)

$(OBJ)/%.o: $(SRC)/*.f $(SRC)/*.F90 $(SRC)/*.c 
	+$(MAKE) -C $(SRC)/

	rm -f $(BASE)/output/tmp/*
	@echo '**** IMPORTANT! ****'
	@echo '** Before start the simulations:'
	@echo '** If in HPC? do not forget to run> ' $(MSG)

clean:
	rm -f $(OBJ)/*.o
	rm -f $(BASE)/output/tmp/*.*

#RS_rt3.o rt3.o radtran3.o radutil3.o radscat3.o radintg3.o radmat.o scat_utilities.o dsd_utilities.o Fresnel_surf.o RT_ncdf_IO.o

