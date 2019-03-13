CF = gfortran

BASE = /home/pga082/GFI/scripts/RT3_NCDF
OBJ = $(BASE)/obj
BIN = $(BASE)/bin
SRC = $(BASE)/src

$(BIN)/RS_RT3_run: $(OBJ)/*.o
	$(CF) -fopenmp -o $@ $(OBJ)/RS_rt3.o $(OBJ)/rt3.o $(OBJ)/radtran3.o \
        $(OBJ)/radutil3.o  $(OBJ)/radscat3.o $(OBJ)/radintg3.o $(OBJ)/radmat.o \
	$(OBJ)/scat_utilities.o $(OBJ)/dsd_utilities.o $(OBJ)/Fresnel_surf.o \
	$(OBJ)/RT_ncdf_IO.o $(OBJ)/mitt_time.o -L/usr/lib -lnetcdff

$(OBJ)/*.o: $(SRC)/*.f $(SRC)/*.F90 $(SRC)/*.c 
	+$(MAKE) -C $(OBJ)/

	rm -f $(BASE)/output/tmp/*

clean:
	rm -f $(OBJ)/*.o
	rm -f $(BASE)/output/tmp/*.*

#RS_rt3.o rt3.o radtran3.o radutil3.o radscat3.o radintg3.o radmat.o scat_utilities.o dsd_utilities.o Fresnel_surf.o RT_ncdf_IO.o

