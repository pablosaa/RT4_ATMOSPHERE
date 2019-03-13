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


---
For more info contact:

Pablo Saavedra G. (pablo.saa@uib.no)

Geophysical Institute, University of Bergen.

(c) 2018 Pablo Saavedra G.

See LICENSE.TXT
