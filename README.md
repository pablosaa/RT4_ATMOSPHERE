# RT4_ATMOSPHERE

A Microwave Radiative Transfer for the Atmosphere

This repository is based on the Radiative Transfer code RT3 and RT4 by Evans et al. 2004?

The input for the code constructed from radiosone data obtained from the Wyoming University and processed by the repository https://github.com/pablosaa/WyoSondes

Then a netCDF file is generate as input for the radiative transfer code.

The code support multiple atmospheric profile entries in a (lat,lon) grid and for a range of time series. It also support multiple microwave frequencies in a single run.

Input parameters are specified by means of a name-list sumarized in the ``input`` file, if that file is not provided then the parameters need to be introduced manually. An example input file is provided at the directory ``input_exemplar/input_base`` containing parameters for simulations suited for the HATPRO radiometer as example. 

The output of the simulations are stored at the directory ``output/TB`` as netCDF files too. The output netCDF contains the following variables:
* Brightness Temperatures for all frequencies and H- & V-polarization at the TOA and ground level,
* The original atmospheric profiles,
* The microphysics profiles, e.g. atmospheric attenuation by gases, clouds, rain
* Surface meteorological variables.

## Compilation
To compile the code run from the base directory: 
  
    > make
  
this will generate the executable file at ``./bin/RS_RT3_run``

## Executing
To run simulations simple execute the generated executable from the ``bin`` directory as:

    {base_dir}/bin> ./RS_RT3_run
  
once the simulations are finished the output file is located at ``output/TB`` and has the following name:
    
    {base_dir}/output/TB/RT3TB_[microphysics_parameters]_fxxxx.nc
  
where ``[microphysics_parameters]`` is a sequence of specifications regarding the parametrization for atmospheric hydrometeors like rain, graupel, snow and ice. The sufix ``_fxxxx`` is just an indicator that the file contains simulation for multiple frequencies.

# Structure of directories
This repository is distributed as follows:

    ./
    ├── bin
    ├── input_exemplar
    ├── obj
    ├── output
    │   ├── TB
    │   └── tmp
    ├── scripts
    └── src

8 directories where ``./bin`` contains the binary executables and input parameter file, ``./src`` has all the source codes, ``./obj`` host the object files from the source codes, ``./output/TB`` is where netCDF output files with the simulations is stored, ``./output/tmp`` is a temporal folder to host auxiliary files during running time, and ``./input_exemplar`` can contain examples of input parameter file with different set of configurations. 

# Data structure
## Dimensions
The database stored in a netCDF file has the following dimensions: time, theta_z, stokes, level, x_grid, y_grid

``theta_z`` is the zenithal angle with 0 deg for Zenith and 180 deg for Nadir. The brightess temperatures angular dimension corresponds to the ``theta_z`` values for the downwelling radiation (_DN_) and ``180 - theta_z`` values for the upwelling radiation (_UP_).

---
For more info contact: pablo.saa@uib.no<br>

(c) 2018, Pablo Saavedra Garfias<br>
Geophysical Institute,<br> University of Bergen,<br> Norway.<br>
See LICENSE
