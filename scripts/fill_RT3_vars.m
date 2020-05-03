% script to fill the netCDF RT3 output files with missed variables:

if exist('OCTAVE_VERSION','builtin'),
  pkg load netcdf;
end

faktor = [1, 1e-2, 100];  % units converting factor, [K, hPa, %]

varsout = {'T2m','P2m','RH2m'};
varsin = {'air_temperature_2m','surface_air_pressure','relative_humidity_2m'};

infile = '/home/pablo/GFI/data/AROME/arome_arctic_full_2_5km_20191231T21Z.nc';
outfile = '/home/pablo/GFI/Repos/RT4_ATMOSPHERE/output/TB/RT3TB_x100-102y170-172_arome_arctic_full_2_5km_20191231T21Z.nc';

if ~exist(infile, 'file'),
  error(['Not found: ' infile]);
end

info = ncinfo(infile);
NTIME = info.Dimensions(1).Length;
for i=1:length(varsin),
  VAR = ncread(infile, varsin{i}, [100, 170, 1, 1], [3, 3, 1, NTIME]);
  ncwrite(outfile, varsout{i}, faktor(i)*squeeze(VAR) );
end

% end of script
