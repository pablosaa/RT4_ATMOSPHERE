CF = gfortran
CC = gcc

LFLAG = -c -fno-automatic

BASE = /home/pga082/GFI/scripts/RT3_NCDF
OBJ = $(BASE)/obj
BIN = $(BASE)/bin
SRC = $(BASE)/src


all: RS_rt3.o rt3.o radtran3.o radutil3.o radscat3.o radintg3.o radmat.o scat_utilities.o dsd_utilities.o Fresnel_surf.o RT_ncdf_IO.o mitt_time.o
	echo "Subroutines compiled..."
#

##RS_RT3_run: RS_rt3.o rt3.o radtran3.o radutil3.o radscat3.o radintg3.o radmat.o scat_utilities.o dsd_utilities.o Fresnel_surf.o RT_ncdf_IO.o
#
##	gfortran -fopenmp -o RS_RT3_run RS_rt3.o rt3.o radtran3.o\
        radutil3.o  radscat3.o radintg3.o radmat.o scat_utilities.o\
          dsd_utilities.o Fresnel_surf.o RT_ncdf_IO.o mitt_time.o -L/usr/lib -lnetcdff
#
RS_rt3.o: $(SRC)/RS_rt3.f $(SRC)/parameters.inc
	gfortran -c -fno-automatic $(SRC)/RS_rt3.f
#
rt3.o: $(SRC)/rt3.f
	gfortran -c -fno-automatic $(SRC)/rt3.f
#
Fresnel_surf.o: $(SRC)/Fresnel_surf.f
	gfortran -c -fno-automatic $(SRC)/Fresnel_surf.f
#
radtran3.o: $(SRC)/radtran3.f
	gfortran -c -fno-automatic $(SRC)/radtran3.f
#
radutil3.o: $(SRC)/radutil3.f
	gfortran -c -fno-automatic $(SRC)/radutil3.f
#
radscat3.o: $(SRC)/radscat3.f
	gfortran -c -fno-automatic $(SRC)/radscat3.f
#
radintg3.o: $(SRC)/radintg3.f
	gfortran -c -fno-automatic $(SRC)/radintg3.f
#
radmat.o: $(SRC)/radmat.f
	gfortran -c -fno-automatic $(SRC)/radmat.f
#
scat_utilities.o: $(SRC)/scat_utilities.f
	gfortran -c -fno-automatic $(SRC)/scat_utilities.f
#
dsd_utilities.o: $(SRC)/dsd_utilities.f
	gfortran -c   $(SRC)/dsd_utilities.f
#
mitt_time.o: $(SRC)/mitt_time.c
	gcc -c $(SRC)/mitt_time.c
#
RT_ncdf_IO.o: $(SRC)/RT_ncdf_IO.F90
	gfortran -c $(SRC)/RT_ncdf_IO.F90 -I/usr/include

