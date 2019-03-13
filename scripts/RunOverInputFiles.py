#! /usr/bin/python3

# script to list all Radiosonde netCDF files in a given directory
# to create the name-list variables for input to the RT3 code.
# This script only change the string for input_file with netCDF files
# located at the specified directory, whereas the other parameters are
# copied either from any specified exemplar input file or by default
# from the existing '../input_exemplar/input_base' sample.
# The script can scan netCDF files recursively with the option -R at
# the end of the command line. See USAGE
#
# (c) 2018 Pablo Saavedra G.
# Geophysical Institute, University of Bergen
# SEE LICENSE.TXT
# ----------------------------------------------------------------- 
#
# USAGE:
# >./RunOverInputFiles arg1 [arg2] [-R]
# 
# the [] means optional parameters,
# Where:
# * arg1: the full Path to directory where the netCDF input files are located,
# * arg2: an ASCII input file containing the parameter's namelist as input
# * -R  : if specified then the scan for netCDF files fill be recursive from arg1
#
# If not arg2 is specified, then the script uses a standard input model for 
# parameters suitable for RPG HATPRO radiometer. 

import glob, os, sys

# Default values:
PATH_BIN = '../bin/'
CMD = PATH_BIN+'RS_RT3_run'
REC = ''
ext = 'nc'
inputpar = '../input_exemplar/input_base'
outputpar = PATH_BIN+'input'
pathin = './'
pathout = ''

if '-R' in sys.argv:
	print('Warning: will do it recursively!')

if len(sys.argv)==4:
	pathin = str(sys.argv[1])
	inputpar = str(sys.argv[2])
	REC = str(sys.argv[3])
elif len(sys.argv)==3:
	inputpar = str(sys.argv[2])
	pathin = str(sys.argv[1])
elif len(sys.argv)==2:
	pathin = str(sys.argv[1])
else:
	print('USAGE:	')
	print('> ./RunOverInputFiles /Disk/path_to_netCDF_files input_base -R')
	print('> ./RunOverInputFiles /Disk/path_to_netCDF_files input_base')
	print('> ./RunOverInputFiles /Disk/path_to_netCDF_files')

if 'R' in REC:
	REC = '**/'

for infile in glob.iglob(pathin+REC+'*.nc'):
	print(infile)
	with open(outputpar,"w") as out:
		with  open(inputpar) as tmp:
			for line in tmp:
				if 'input_file=' in line:
					out.write("	input_file='"+infile+"',\n")
				else:
					out.write(line)

	os.system('head '+outputpar) #CMD)
	os.system('rm -f '+outputpar)
print('End RT over all netCDF files. ')
