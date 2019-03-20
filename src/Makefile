CF = gfortran
CC = gcc

LFLAG = -c -fno-automatic -std=legacy
OFLAG = -o $@

#BASE = .. #/home/pga082/GFI/scripts/RT3_NCDF
OBJ = $(BASE)/obj
BIN = $(BASE)/bin
SRC = $(BASE)/src

#RS_rt3.o rt3.o radtran3.o radutil3.o radscat3.o radintg3.o radmat.o scat_utilities.o dsd_utilities.o Fresnel_surf.o RT_ncdf_IO.o mitt_time.o

#all: $(OBJ)/*.o 
#

$(BIN)/RS_RT3_run: $(OBJ)/RS_rt3.o $(OBJ)/rt3.o $(OBJ)/radtran3.o $(OBJ)/radutil3.o $(OBJ)/radscat3.o $(OBJ)/radintg3.o $(OBJ)/radmat.o $(OBJ)/scat_utilities.o $(OBJ)/dsd_utilities.o $(OBJ)/Fresnel_surf.o $(OBJ)/RT_ncdf_IO.o
	echo "Subroutines compiled..."

##	gfortran -fopenmp -o RS_RT3_run RS_rt3.o rt3.o radtran3.o\
        radutil3.o  radscat3.o radintg3.o radmat.o scat_utilities.o\
          dsd_utilities.o Fresnel_surf.o RT_ncdf_IO.o mitt_time.o -L/usr/lib -lnetcdff
#
$(OBJ)/RS_rt3.o: $(SRC)/RS_rt3.f $(SRC)/parameters.inc
	gfortran $(LFLAG) $(OFLAG) $(SRC)/RS_rt3.f
#
$(OBJ)/rt3.o: $(SRC)/rt3.f
	gfortran $(LFLAG) $(OFLAG) $(SRC)/rt3.f
#
$(OBJ)/Fresnel_surf.o: $(SRC)/Fresnel_surf.f
	gfortran $(LFLAG) $(OFLAG) $(SRC)/Fresnel_surf.f
#
$(OBJ)/radtran3.o: $(SRC)/radtran3.f
	gfortran $(LFLAG) $(OFLAG) $(SRC)/radtran3.f
#
$(OBJ)/radutil3.o: $(SRC)/radutil3.f
	gfortran $(LFLAG) $(OFLAG) $(SRC)/radutil3.f
#
$(OBJ)/radscat3.o: $(SRC)/radscat3.f
	gfortran $(LFLAG) $(OFLAG) $(SRC)/radscat3.f
#
$(OBJ)/radintg3.o: $(SRC)/radintg3.f
	gfortran $(LFLAG) $(OFLAG) $(SRC)/radintg3.f
#
$(OBJ)/radmat.o: $(SRC)/radmat.f
	gfortran $(LFLAG) $(OFLAG) $(SRC)/radmat.f
#
$(OBJ)/scat_utilities.o: $(SRC)/scat_utilities.f
	gfortran $(LFLAG) $(OFLAG) $(SRC)/scat_utilities.f
#
$(OBJ)/dsd_utilities.o: $(SRC)/dsd_utilities.f
	gfortran $(LFLAG) $(OFLAG) $(SRC)/dsd_utilities.f
#
$(OBJ)/mitt_time.o: $(SRC)/mitt_time.c
	gcc -c $(OFLAG) $(SRC)/mitt_time.c
#
$(OBJ)/RT_ncdf_IO.o: $(SRC)/RT_ncdf_IO.F90
	gfortran -c $(OFLAG) $(SRC)/RT_ncdf_IO.F90 -I$(INCPATH)
