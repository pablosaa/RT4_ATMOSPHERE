! *****************************************************************************************
! Set of subroutines to manage netCDF input/output files.
! -----------------------------------------------------------------------------------------
! * module nctoys: Contains axiliary modules for managing netCDF files. Interfaces to
!                  external functions used by netCDF subroutines.
! ---
! * read_arome   : read netCDF files from AROME-Arctic output files,
! * read_wrf     : read netCDF files from WRF output files,
! * read_wyosonde: read netCDF files from radiosondes with homoginized layers,
! * createncdf   : create and setup variables and dimensions for RT output files,
! * storencdf    : save the TB variables into the previous created netCDF file,
! * MP_storencdf : save the microphysics variables into the created netCDF file.
!
! (c) 2018 Pablo Saavedra G. (pablo.saa@uib.no)
! Geophysical Institute, University of Bergen
! SEE LICENSE.TXT
! --------------------------------------------------------------------------------------

module nctoys
  public
  type ncname
     character(len=15), allocatable, dimension(:) :: vars
     character(len=15), allocatable, dimension(:) :: dims
  end type ncname

  ! Interface for generic function for UnixTime and Date conversions:
  ! How to use (3 following options):
  ! > unixtime = getUnixTime(date)
  ! > unixtime = getUnixTime(DateString)
  ! > date(i)  = getUnixTime(unixtime(i) )
  !
  ! Where:
  ! * date: is a N x 6 array with rows as [year, month, day, hour, min, sec]
  ! * DateString: is a N x 19 char with row as "YYYY-MM-DD_hh:mi:se"
  ! * unixtime  : a 1D array with elements as days since epoch.
  ! * i         : indice to indiates one element at a time
  interface getUnixTime
     module procedure getF2UnixTime, getCH2UnixTime, getUnixTime2Date
  end interface getUnixTime
  ! -- end of interface UnixTime generic function
  
contains


  ! ----
  ! Subroutine to get dimensions from netCDF file needed for RT
  ! Four dimension are mandatory in the following order:
  ! [gridx, gridy, layers, time]
  !
  subroutine get_dims_allocate_vars(ncid, names)
    use netcdf
    use variables, only : nx_in, nx_fin, ny_in, ny_fin, ngridx, ngridy, nlyr, ntime
    implicit none
    integer, intent(in) :: ncid
    type(ncname), intent(in) :: names
    ! local variables:
    integer :: i, status, dim_id, dim_len
    
    do i = 1, size(names%dims)
       status = nf90_inq_dimid(ncid, trim(names%dims(i)), dim_id)
       call check_nc(status, 'dimension name does not exist!', .TRUE.)
       status = nf90_inquire_dimension(ncid, dim_id, len = dim_len )     
       call check_nc(status, 'not possible to retrieve dimensions', .TRUE.)

       select case(i)
       case(1)
          if(nx_in.EQ.0) nx_in = 1
          if(nx_fin.EQ.0) nx_fin = dim_len
       case(2)
          if(ny_in.EQ.0) ny_in = 1
          if(ny_fin.EQ.0) ny_fin = dim_len
       case(3)
          nlyr = dim_len
       case(4)
          ntime = dim_len
       case default
          print*, 'netCDF dimension '//trim(names%dims(i))//' unused!'
       end select
    end do
    ngridx = nx_fin - nx_in +1
    ngridy = ny_fin - ny_in +1

    ! Allocate global variables according to dimensions:
    CALL allocate_RT3_variables

  end subroutine get_dims_allocate_vars

  ! ----
  ! Subroutine to check status of netCDF operations:
  subroutine check_nc(status, message, FLAG)
    use netcdf
    implicit none
    integer, intent(in) :: status
    character(len=*), optional, intent(in) :: message
    logical, optional, intent(in) :: FLAG
    if(status.EQ.NF90_NOERR) return
    ! Ups... get this point cause error happened?
    print *, NF90_STRERROR(status)
  
    if(present(message) ) print*, message
    if(present(FLAG) ) stop
    return
  end subroutine check_nc

  ! ----
  ! Subroutine to work with converting Date to unix-time (uses external functions)
  ! ---
  function getUnixTime2Date(unixtime) result(datum)
    use, intrinsic :: iso_c_binding
    implicit none
    interface
       subroutine UnixTime2Date(unixt, datum) bind(c, name='UnixTime2Date')
         import
         real(kind=c_double), value :: unixt
         integer(kind=c_int) :: datum(6)
       end subroutine UnixTime2Date
    end interface
    real(kind=8) :: unixtime
    integer :: datum(6)

    call UnixTime2Date(unixtime, datum)
  end function getUnixTime2Date

  ! ----
  function getF2UnixTime(datum) result(unixtime)
    use, intrinsic :: iso_c_binding

    implicit none
    ! Interface to the C code for Unix time retrieval:
    interface
       subroutine F2UnixTime(ntime, datum, val) bind(c, name='F2UnixTime')
         import 
         integer(kind=c_int), value :: ntime
         integer(kind=c_int) :: datum(ntime, 6)
         real(kind=c_double) :: val(ntime)
       end subroutine F2UnixTime
    end interface
    
    integer(c_int), intent(in) :: datum(:,:)
    real(c_double), allocatable :: unixtime(:)

    integer :: ntime
    ntime = size(datum,1)
    allocate(unixtime(ntime) )
    call F2UnixTime(ntime, datum, unixtime)
    return
    
  end function GetF2UnixTime

  function getCH2UnixTime(datum) result(unixtime)
    use, intrinsic :: iso_c_binding

    implicit none
    ! Interface to the C code for Unix time retrieval:
    interface
       subroutine CH2UnixTime(ntime, nlen,  datum, val) bind(c, name='CH2UnixTime')
         import 
         integer(kind=c_int), value :: ntime, nlen
         character(kind=c_char,len=1) :: datum(ntime)
         real(kind=c_double) :: val(ntime)
       end subroutine CH2UnixTime
    end interface
    character(kind=c_char, len=*), intent(in) :: datum(:)
    real(c_double), allocatable :: unixtime(:)

    integer :: ntime, nlen
    nlen = len(datum)
    ntime = size(datum,1)
    allocate(unixtime(ntime) )
    call CH2UnixTime(ntime, nlen, datum, unixtime)
    return
    
  end function GetCH2UnixTime

end module nctoys

module meteo_tools
  implicit none
contains
  ! --------------------------------------------------------------------
  ! Subroutine to interpolate extra observations angles in case the
  ! input value is out of limits, extrapolation is attended with the
  ! two first (lower extrapolation) or last (upper extrapolation)
  ! values of the vector.
  ! ---
  subroutine interp1_1D(Xin, Yin, Nin, X0, N0, Xout, Nout, Yout, Nstk, Nflx, Nlev)
  
  ! Interpolated abssice values must increase monotonically, 
  ! the ordenate values can be decreasing or increasing
  ! Repeated values in Xin and X0 are omitted
  
  implicit none
  integer, intent(in) :: Nin, N0, Nstk, Nflx, Nlev, Nout
  real(kind=8), intent(in) :: Xin(Nin), X0(N0), Yin(Nin,Nstk,Nflx,Nlev)
  real(kind=8), intent(out):: Xout(Nout)
  real(kind=8), intent(out):: Yout(Nout,Nstk,Nflx,Nlev)
  integer :: h, idx, idn, iout
  real(kind=8) :: a, b, Fa(Nstk, Nflx, Nlev), Fb(Nstk, Nflx, Nlev), F0(Nstk, Nflx, Nlev)
  real(kind=8) :: eps = 1.0E-4
  
  ! union of Xin and X0 sets
  Yout = -69.8
  Yout(1:Nin,:,:,:) = Yin
  
  Xout = -69.0
  Xout(1:Nin) = Xin
    
  do h = 1, N0
     !if(any(abs(Xout-X0(h)).LT.eps)) cycle
     !iout = iout + 1
     if(X0(h).GT.maxval(Xout)) then
        idx = maxloc(Xout,1)+1
        a = Xout(idx-2)
        b = Xout(idx-1)
        Fa = Yout(idx-2,:,:,:)
        Fb = Yout(idx-1,:,:,:)
     else if(X0(h).LT.minval(Xout)) then
        idx = 1
        a = Xout(idx)
        b = Xout(idx+1)
        Fa = Yout(idx,:,:,:)
        Fb = Yout(idx+1,:,:,:)
     else
        idx = minloc(Xout,DIM=1,MASK=X0(h)<=Xout)
        a = Xout(idx-1)
        b = Xout(idx)
        Fa = Yout(idx-1,:,:,:)
        Fb = Yout(idx,:,:,:)
     end if
     
     if(abs(X0(h)-b).LT.eps) then
        F0 = Fb
     else if(abs(X0(h)-a).LT.eps) then
        F0 = Fa
     else
        F0 = Fa + (X0(h)-a)*(Fb-Fa)/(b-a)
     end if
     idn = Nin+h-1
     Xout(idx+1:idn+1) = Xout(idx:idn)
     Xout(idx) = X0(h)

     Yout(idx+1:idn+1, :, :, :) = Yout(idx:idn,:,:,:)
     Yout(idx, :, :, :)   = F0
     
  end do
  
end subroutine interp1_1D
! ----/

! ---------------------------------------------------------------
! ELEMENTAL FUNCTION Convert mixing ratio to Relative Humidity.
! HUMIDITY CONVERSION FORMULAS
! Calculation formulas for humidity (B210973EN-F)
! By (c) VAISALA 2003
! https://www.hatchability.com/Vaisala.pdf
!
! ---
! (c) 2019 Pablo Saavedra G.
! Geophysical Institute, University of Bergen
! See LICENSE
!
! ---
! -> MIXR: Vapour mixing ratio [kg/kg]
! -> P   : Pressure [hPa]
! -> T   : Temperature [K]
! <- RH  : Relative Humidity [%]
elemental function mixr_to_rh(MIXR, P, T) result(RH)
  implicit none
  real(kind=8), intent(in) :: MIXR, P, T
  real(kind=8) :: RH

  integer :: i, j
  real :: PWS, etha, A
  real, parameter :: COEFF = 2.16679 ! [g K J^-1]
  real, parameter :: Tc = 647.096  ! critical temperature [K]
  real, parameter :: Pc = 220640   ! critical pressure [hPa]
  real, parameter :: B  = 0.6219907 ! constant for air [kg/kg]
  real, parameter :: CC(6) = (/ -7.85951783, 1.84408259, &
       & -11.7866497, 22.6807411, -15.9618719, 1.80122502 /)
  real, parameter :: EE(6) = (/1.0, 1.5, 3.0, 3.5, 4.0, 7.5/)

  etha = 1.0d0 - T/Tc
  A = 0.0d0
  do i=1, 6
     A = A + CC(i)*(etha**EE(i))
  end do
  PWS = Pc*exp(A*Tc/T)
  !RH = 1E-3*MIXR*T/PWS/COEFF
  !RH = 100*PW/PWS
  RH = 100*MIXR*P/(MIXR + B)/PWS
  return
end function mixr_to_rh
! ----/

! -------------------------------------------------------------------
! ELEMENTAL Function to convert Specific [kg/kg] to
! mixing rations [kg/kg] quantitis.
! -> Q_x  : Specific content in [kg/kg]
! <- mixr : Mixing ratio in [kg/kg]
! ---
elemental function qx_to_mixr(Q_x) result(MIXR)
  implicit none
  real(kind=8), intent(in) :: Q_x
  real(kind=8) :: MIXR

  MIXR = Q_x/(1.0-Q_x)
  return
end function qx_to_mixr
! ----/

! _________________________________________________________________
! ELEMENTAL FuNCTION Specific Humidity to Relative Humidity
! -> Qv : Specific Humidity [kg/kg]
! -> P  : Pressure [hPa]
! -> T  : Temperature [K]
! <- RH : Relative Humidity [%]
! ---
elemental function qv_to_rh(QV, P, T) result(RH)
  implicit none
  real(kind=8), intent(in) :: QV, P, T
  real(kind=8) :: RH

  ! Local variables:
  real(kind=8) :: MIXR
  
  ! Converting Specific Humidity to Mixing Ratio:
  MIXR = qx_to_mixr(QV)
  
  RH = mixr_to_rh(MIXR, P, T)

  return
end function qv_to_rh
! ----/

! ____________________________________________________________________
! ELEMENTAL FUNCTION TO CONVERT WRF perturbation potential temperature
! (theta-t0) to Temperature.
!
! Based on "[Wrf-users] perturbation potential temperature"
! https://mailman.ucar.edu/pipermail/wrf-users/2010/001896.html
! Where thete: Potential Temperature and T0=300K
! kappa: the Poisson constant (kappa = R/c_p), the ratio of the gas
! constant R to the specific heat at constant pressure c_p.
!
! -> Tper : Perturbation potential temperature [K]
! -> P    : Pressure [hPa]
! <- T    : Temperature [K]
!
! ---
elemental function PERTHETA2T(Tper, P) result(T)
  implicit none
!  integer, intent(in) :: nx, ny, nz, nt
  real(kind=8), intent(in) :: Tper, P
  real(kind=8) :: T

  real(kind=8) :: theta
  real, parameter :: T0 = 300.0  ! [K]
  real, parameter :: P0 = 1000.0 ! standard pressure [hPa]
  real, parameter :: kappa = 0.2854 !  For dry air

  theta = Tper + T0    ! Potential Temperature [K]
  T = theta*(P/P0)**kappa   ! [K]
  return
end function PERTHETA2T
! ----/

! __________________________________________________________________
! Subroutine to convert Wind U and V component into
! wind speed, and direction
! -> U, V : Wind components U and V [m/s]
! <- WS   : Wind speed [m/s]
! <- WD   : Wind direction [deg]
!
! (c) 2019, Pablo Saavedra G.
! Geophysical Institute, University of Bonn
! ---
elemental subroutine Wind_UV2speeddir(U, V, WS, WD)
  implicit none
  real(kind=8), intent(in) :: U, V
  real(kind=8), intent(out) :: WS, WD
  real, parameter :: PI2deg = 45.0/atan(1.0)
  
  WS = sqrt( U*U + V*V)
  WD = modulo(360.0 - atan2(U, V)*PI2deg, 360.0)
  return
end subroutine Wind_UV2speeddir
! ----/


! _______________________________________________________________________
! Subroutine to calculate the Geopotential Height given T, P, Qv
!
! - INPUT:
! * T(NX, NY, NZ, NT)    : Temperature [K]
! * P(NX, NY, 0:NZ, NT)  : Pressure [hPa]
! * Qv(NX, NY, NZ, NT)   : Specific Humidity [kg/kg]
! - OUTPUT:
! * Z(NX, NY, NZ, NT)    : Geopotential Height [km] 
!
! ---
! (c) 2020, Pablo Saavedra G.
! Geophysical Institute, University of Bergen
! See LICENSE
! -----------------------------------------------------------------------
function Calculate_GeopotentialZ(T, P, Qv) result(Z) 
  implicit none
  real(kind=8), intent(in), dimension(:,:,:,:) :: T, P, Qv
  real(kind=8), allocatable :: Z(:,:,:,:)

  ! ** local variables
  integer :: i, NZ
  real(kind=8), allocatable, dimension(:,:,:,:) :: MIXR, TV, del_P, TVdP !
  real, parameter :: Rd = 287.0       ! [J/kg/K] dry air gass constant
  real, parameter :: g0 = 9.807       ! [m/s^2]  gravity constant

  MIXR = qx_to_mixr(Qv)   ! [kg/kg]
  TV = VirtualTemperature(T, MIXR)

  NZ = size(T, 3)
  if(size(P,3).NE.NZ+1) stop 'Layer dimension for P needs to be 1 nore than Z'
  del_P = P(:,:,1:NZ,:) - P(:,:,2:NZ+1,:)

  ! temporal integrable variable = Tv/P *dP
  TVdP = TV*del_P/P(:,:,1:NZ,:)

  ! Initializing Z = 0. source T
  Z = T*0.0d0
  do i = 1, NZ
     Z(:,:,i,:) = sum(TVdP(:,:,1:i,:), dim=3)
  end do

  Z = (1.0E-3)*Z*Rd/g0  ! converting to [km]

  return
end function Calculate_GeopotentialZ
! -----/

! _______________________________________________________________________
! Subroutine to calculate the Virtual Temperature given T, MIXR
!
! --------------------------------------------------------------------
! Function convert Temperature to Virtual Temp.
! -> T    : Temperature [K]
! -> MIXR : Vapour Mixing ratio [kg/kg]
! <- Tv   : Virtual Temperature [K]
! ---
! (c) 2020, Pablo Saavedra G.
! Geophysical Institute, University of Bergen
! See LICENSE
! ---
elemental function VirtualTemperature(T, MIXR) result (TV)
  implicit none
  real(kind=8), intent(in) :: T, MIXR
  real(kind=8) :: TV
  ! ** local variables
  ! Rd gas constant dry air, Rv gas constant water vapour:
  real, parameter :: Rd = 287, Rv = 461.5  ! [J/kg/K]
  real, parameter :: eps = Rd/Rv           ! ~ 0.622

  TV = T*(1.0 + MIXR/eps)/(1.0 + MIXR)
  return
end function VirtualTemperature
! ----/

! --------------------------------------------------------------------
! Function to convert kg/kg to kg/m^3
! -> Qx  : Specific quantity e.g. humidity [kg/kg]
! -> T   : Temperature [K]
! -> P   : Pressure [hPa]
! -> MIXR: vapour mixing ration [kg/kg]
! <- RHOx: Specific quantity [kg/m^3]
! ---
elemental function MassRatio2MassVolume(Q_x, T, P ,MIXR) result(RHO_x)
  implicit none
  real(kind=8), intent(in) :: Q_x, T, P, MIXR
  real(kind=8) :: RHO_x
  ! Local variables:
  real(kind=8) :: Tv, RHO_air
  real, parameter :: Rd = 287, Rv = 461.5  ! [J/kg/K]
  
  Tv = VirtualTemperature(T, MIXR)
  RHO_air = (1.0E2*P)/Tv/Rd  ! [kg/m^3]
  RHO_x   = Q_x*RHO_air      ! [kg/m^3]
  
  return
end function MassRatio2MassVolume
! ----/

! ---------------------------------------------------------------------
! Converting hydrometeor variable units from kg/kg to g/m^3
!
subroutine Convert_Hydrometeors_units(Is_Specific_VAR)
  !use meteo_tools, only: MassRatio2MassVolume, qx_to_mixr
  use variables, only: temp_tmp, press_tmp, cloud_water_tmp, &
       rain_water_tmp, cloud_ice_tmp, &
       snow_tmp, graupel_tmp, mixr_tmp
  implicit none
  logical, intent(in) :: Is_Specific_VAR

  if(Is_Specific_VAR) then
     ! TRUE if variables are specific quantities and 
     ! need to be converted to mixing rations:
     cloud_water_tmp = qx_to_mixr(cloud_water_tmp)
     rain_water_tmp  = qx_to_mixr(rain_water_tmp)
     cloud_ice_tmp   = qx_to_mixr(cloud_ice_tmp)
     snow_tmp = qx_to_mixr(snow_tmp)
     graupel_tmp = qx_to_mixr(graupel_tmp)
  end if
  
  cloud_water_tmp = 1.0E3*MassRatio2MassVolume(cloud_water_tmp, &
       & temp_tmp(:,:,1:,:), press_tmp(:,:,1:,:), mixr_tmp )

  rain_water_tmp = 1.0E3*MassRatio2MassVolume(rain_water_tmp, &
       temp_tmp(:,:,1:,:), press_tmp(:,:,1:,:), mixr_tmp )

  cloud_ice_tmp = 1.0E3*MassRatio2MassVolume(cloud_ice_tmp, &
       temp_tmp(:,:,1:,:), press_tmp(:,:,1:,:), mixr_tmp )

  snow_tmp = 1.0E3*MassRatio2MassVolume(snow_tmp, &
       temp_tmp(:,:,1:,:), press_tmp(:,:,1:,:), mixr_tmp )

  graupel_tmp = 1.0E3*MassRatio2MassVolume(graupel_tmp, &
       temp_tmp(:,:,1:,:), press_tmp(:,:,1:,:), mixr_tmp )

end subroutine Convert_Hydrometeors_units
! -----/

end module meteo_tools


! ______________________________________________________________________________________
! --------------------------------------------------------------------------------------
! Subroutine to read AROME-Arctiv simulation outputs from NetCDF-file generated by
! met.no ARMORE-Arctic public data in their web archive.
!
! The loaded variables are introduced to the RT3/RT4 radiative transfer code for
! calculation of Brightness Temperature Fields.
! -------------------------------------------------------------------------------------
subroutine read_arome(ncflen, ncfile, del_xy, origin_str)
  use netcdf
  use variables
  use nctoys
  use meteo_tools
  
  implicit none

  integer, intent(in) :: ncflen
  character(len=ncflen), intent(in) :: ncfile
  ! -- Here come variables declaration alike the onses used by RT3
  real(kind=4), intent(out) :: del_xy(2)
  character(len=*), intent(out) :: origin_str

  ! -- end of variable declaration
  integer, parameter :: NVarIn = 18, Nindim = 4
  integer :: status
  integer :: ncid, ndims_in, nvars_in ,ngatts_in, unlimdimid_in
  integer :: VarID, dim_len

  character(len=30) :: varname
  character(len=10) :: dim_name
  ! Auxiliary variables for magnitude conversions:
  integer :: i, N_tempX, N_tempY
  real :: faktor
  real(kind=8), allocatable, dimension(:,:,:,:) :: U_Vel, V_vel, QV
  real, allocatable, dimension(:) :: sigma_hybrid
  
  character(len=69), dimension(*), parameter :: arome_name = & !NVarIn
       [ character(len=69):: &
       'surface_geopotential', &
       'air_temperature_2m', &
       'surface_air_pressure', &
       'specific_humidity_2m', &
       'latitude', &
       'longitude', &
       'atmosphere_boundary_layer_thickness', &
       'hybrid', &
       'air_temperature_ml', &
       'specific_humidity_ml', &
       'mass_fraction_of_cloud_condensed_water_in_air_ml', &
       'mass_fraction_of_rain_in_air_ml', &
       'mass_fraction_of_cloud_ice_in_air_ml', &
       'mass_fraction_of_snow_in_air_ml', &
       'mass_fraction_of_graupel_in_air_ml', &
       'x_wind_ml', &
       'y_wind_ml', &
       'time']

  type(ncname) :: aromename

  ! -----------------------------------------------------------------------------
  ! Open file and read directory
  print*,'AROME-Arctic netCDF input files is : ', trim(ncfile)
  
  status = nf90_open(ncfile,NF90_NOWRITE,ncid)
  call check_nc(status, 'Cannot open AROME netCDF-file '//trim(ncfile), .TRUE.)


  aromename = ncname(dims=[ character(len=6):: &
       & 'x', &
       & 'y', &
       & 'hybrid', &
       & 'time'])

  ! Obtain dimensions and allocate RT variables:
  CALL get_dims_allocate_vars(ncid, aromename)


  ! Allocation For local/temporal variables:
  allocate( sigma_hybrid(nlyr) )
  allocate( QV(ngridx, ngridy, 0:nlyr, ntime) )
  allocate( U_Vel(ngridx, ngridy, nlyr, ntime) )
  allocate( V_Vel(ngridx, ngridy, nlyr, ntime) )
  
  do i = 1, size(arome_name)
     print*, 'Reading variable: ', trim(arome_name(i) )
     status = nf90_inq_varid(ncid, trim(arome_name(i) ), VarId)
     call check_nc(status, 'WARNING: for '//trim(arome_name(i)) )

     status = nf90_get_att(ncid, VarId, 'scale_factor', faktor)
     if(status.NE.NF90_NOERR) faktor = 1.

     select case(trim(arome_name(i) ))
        ! *** Reading for AROME SURFACE variables:
     case('surface_geopotential')
        status = nf90_get_var(ncid, VarId, hgt_tmp(:,:,0:0,:), &
             start=(/nx_in, ny_in, 1, 1/), &
             count=(/ngridx, ngridy, 1, ntime/) )
        hgt_tmp(:,:,0,:) = max(1.0E-3*hgt_tmp(:,:,0,:)/9.807, 0.0)
        
     case('air_temperature_2m')
        status = nf90_get_var(ncid, VarId, temp_tmp(:, :, 0:0, :), &
             & start=(/nx_in, ny_in , 1, 1/), count=(/ngridx, ngridy, 1, ntime/) )

     case('surface_air_pressure')
        status = nf90_get_var(ncid, VarId, press_tmp(:, :, 0:0, :), &
             & start=(/nx_in, ny_in , 1, 1/), count=(/ngridx, ngridy, 1, ntime/) )

     case('specific_humidity_2m')
        status = nf90_get_var(ncid, VarId, QV(:, :, 0:0, :), &
             & start=(/nx_in, ny_in , 1, 1/), count=(/ngridx, ngridy, 1, ntime/) )
        
     case('latitude')
        status = nf90_get_var(ncid, VarId, lat, start=(/nx_in, ny_in, 1/), count=(/ngridx, ngridy, 1/))
        
     case('longitude')
        status = nf90_get_var(ncid, VarId, lon, start=(/nx_in, ny_in, 1/), count=(/ngridx, ngridy, 1/))
        
     case('hybrid') ! adapt from sigma to pressure levels
        status = nf90_get_var(ncid, VarId, sigma_hybrid, start=(/1/), count=(/nlyr/))

     case('air_temperature_ml')
        status = nf90_get_var(ncid, VarId, temp_tmp(:, :, nlyr:1:-1, :), &
             & start=(/nx_in, ny_in , 1, 1/), count=(/ngridx, ngridy, nlyr, ntime/))
        temp_tmp(:, :, 1:nlyr, :) = faktor*temp_tmp(:, :, 1:nlyr, :)

     case('specific_humidity_ml')
        status = nf90_get_var(ncid, VarId, QV(:, :, nlyr:1:-1, :), &
             & start=(/nx_in, ny_in , 1, 1/), count=(/ngridx, ngridy, nlyr, ntime/) )
        QV(:, :, 1:nlyr, :) = faktor*QV(:, :, 1:nlyr, :)
        
     case('mass_fraction_of_cloud_condensed_water_in_air_ml')
        status = nf90_get_var(ncid, VarId, cloud_water_tmp(:, :, nlyr:1:-1, :), &
             & start=(/nx_in, ny_in , 1, 1/), count=(/ngridx, ngridy, nlyr, ntime/) )
        
     case('mass_fraction_of_rain_in_air_ml')
        status = nf90_get_var(ncid, VarId, rain_water_tmp(:, :, nlyr:1:-1, :), &
             & start=(/nx_in, ny_in , 1, 1/), count=(/ngridx, ngridy, nlyr, ntime/) )
        rain_water_tmp = faktor*rain_water_tmp
        
     case('mass_fraction_of_cloud_ice_in_air_ml')
        status = nf90_get_var(ncid, VarId, cloud_ice_tmp(:, :, nlyr:1:-1, :), &
             & start=(/nx_in, ny_in , 1, 1/), count=(/ngridx, ngridy, nlyr, ntime/) )
     case('mass_fraction_of_snow_in_air_ml')
        status = nf90_get_var(ncid, VarId, snow_tmp(:, :, nlyr:1:-1, :), &
             & start=(/nx_in, ny_in , 1, 1/), count=(/ngridx, ngridy, nlyr, ntime/) )
     case('mass_fraction_of_graupel_in_air_ml')
        status = nf90_get_var(ncid, VarId, graupel_tmp(:, :, nlyr:1:-1, :), &
             & start=(/nx_in, ny_in , 1, 1/), count=(/ngridx, ngridy, nlyr, ntime/) )
        graupel_tmp = faktor*graupel_tmp
        
     case('x_wind_ml')
        status = nf90_get_var(ncid, VarId, V_Vel(:, :, nlyr:1:-1, :), &
             & start=(/nx_in, ny_in+1 , 1, 1/), count=(/ngridx, ngridy, nlyr, ntime/) )
        V_Vel = faktor*V_Vel
        
     case('y_wind_ml')        
        status = nf90_get_var(ncid, VarId, U_Vel(:, :, nlyr:1:-1, :), &
             & start=(/nx_in+1, ny_in , 1, 1/), count=(/ngridx, ngridy, nlyr, ntime/))
        U_Vel = faktor*U_Vel
        
     case('time')
        status = nf90_get_var(ncid, VarId, UnixTime, start=(/1/), count=(/ntime/) )
        UnixTime = UnixTime/60/60/24
     case('atmosphere_boundary_layer_thickness')
        status = nf90_get_var(ncid, VarId, PBLH, &
             & start=(/nx_in, ny_in , 1, 1/), count=(/ngridx, ngridy, 1, ntime/) )
        
     case default
        print*, 'WARNING: AROME variable ', trim(arome_name(i)),' not recognized.'
     end select
     call check_nc(status, trim(arome_name(i)), .true.)
  end do

  ! 1) Converting sigma_hybrid to pressure levels
  print*, '1. converting sigma to P'
  do i = 1, nlyr
     press_tmp(:, :, nlyr-i+1, :) = sigma_hybrid(i)*press_tmp(:, :, 0, :)
  end do
  press_tmp = press_tmp*1.0E-2  ! [hPa]
  ! --
  ! 2) Converting vapour mixing ratio to Relative Humidity
  print *, '2. converting Qv to RH'
  relhum_tmp = qv_to_rh(QV, press_tmp, temp_tmp)
  
  ! --
  ! 3) Converting specific humidity units. XXto vapour mixing ratio
  print *, '3. converting units for  Qv'
  mixr_tmp = qx_to_mixr(QV(:,:,1:,:))  ! [kg/kg]
  
  ! --
  ! 4) Calculating Geopotential Altitude [km]:
  hgt_tmp(:,:,1:nlyr,:) = Calculate_GeopotentialZ( &
       temp_tmp(:,:,1:nlyr,:), press_tmp, QV(:, :, 1:nlyr,:) )

  ! --
  ! 5) Converting Wind U and V components to Windspeed and Direction:
  print *, '5. converting U V to speed dir'
  call Wind_UV2speeddir(U_Vel, V_Vel, windvel_tmp, winddir_tmp)

  ! --
  ! 6) Converting Specific masses [kg/kg] to mass-volume [g/m^3]:
  call Convert_Hydrometeors_units(Is_Specific_VAR=.true.)
  
  ! --
  qidx = 15
  del_xy = 2.5  ! [km]

  ! ****************************
  ! Retrieving Global Attribute:
  status = nf90_get_att(ncid, NF90_GLOBAL, 'title', origin_str)
  
  ! ****************************
  ! Closing the NetCDF file
  status = nf90_close(ncid)

  deallocate(sigma_hybrid, U_Vel, V_Vel, QV, stat=status)
  if(status/=0) stop 'AROME, something wrong with deallocation?'

  return
end subroutine read_arome
! ---------------------------------------------------------------------/

! _____________________________________________________________________
! ---------------------------------------------------------------------
! Subroutine to read WRF simulation outputs from NetCDF-file generated by standard
! WRF code.
!
! The loaded variables are introduced to the RT3/RT4 radiative transfer code for
! calculation of Brightness Temperature Fields.
! ----------------------------------------------------------------------
subroutine read_wrf(ncflen, ncfile, del_xy, origin_str)
  use netcdf
  use variables
  use nctoys
  use meteo_tools
  
  implicit none

  integer, intent(in) :: ncflen
  character(len=ncflen), intent(in) :: ncfile

  ! -- Here come variables declaration alike the onses used by RT3
  real(kind=4), intent(out) :: del_xy(2)
  character(len=*), intent(out) :: origin_str

  ! -- end of variable declaration

  integer :: status
  integer :: ncid, ndims_in, nvars_in, ngatts_in, unlimdimid_in
  integer :: VarId, dim_len
  character(len=15) :: dim_name
  character(len=30) :: varname
  integer :: i, j, k, NN
  
  real(kind=8), allocatable, dimension(:,:,:,:) :: U_Vel, V_Vel, mixratio, pertur_T
  character(len=:), allocatable, dimension(:) :: TimeStamp
  ! List of variable names to load (WRF convention)
  character(len=15), dimension(17), parameter :: wrfvarname = [ &
       &character(len=15):: &
       & 'T2', &
       & 'PSFC', &
       & 'Q2', &
       & 'XLAT', &
       & 'XLONG', &
       & 'PHB', &
       & 'P_HYD', &
       & 'T', &
       & 'QVAPOR', &
       & 'QCLOUD', &
       & 'QRAIN', &
       & 'QICE', &
       & 'QSNOW', &
       & 'QGRAUP', &
       & 'V', &
       & 'U', &
       & 'Times']

  type(ncname) :: wrfname
  ! --------------------------------------------------------------------
  ! Open file and read directory
  print*,'WRF netCDF input files is', ncflen, ' : ', trim(ncfile)
  
  status = nf90_open(ncfile, NF90_NOWRITE, ncid)
  call check_nc(status, 'Cannot open netCDF-file: '//trim(ncfile), .true.)

  ! Defining dimension names in WRF netCDF file:
  wrfname % dims = [ character(len=11):: &
       & 'west_east', &
       & 'south_north', &
       & 'bottom_top', &
       & 'Time']


  ! Obtain dimensions and allocate RT variables:
  CALL get_dims_allocate_vars(ncid, wrfname)

  ! For local variables:
  allocate(pertur_T(ngridx, ngridy, nlyr, ntime) )
  allocate(mixratio(ngridx, ngridy, 0:nlyr, ntime) )
  allocate(U_Vel(ngridx, ngridy, nlyr, ntime) )
  allocate(V_Vel(ngridx, ngridy, nlyr, ntime) )
  allocate(character(len=19) :: TimeStamp(ntime) )
  
  qidx = 15

  ! loop over all variables needed from WRF netCDF file
  do i=1, size(wrfvarname)
     status = nf90_inq_varid(ncid, trim(wrfvarname(i) ), VarId)
     call check_nc(status, 'Variable '//trim(wrfvarname(i))//' cannot be read.')

     select case(trim(wrfvarname(i) ))
        ! *** Reading for WRF SURFACE variables:
     case('T2')
        status = nf90_get_var(ncid, VarId, temp_tmp(:, :, 0:0, :), &
             start=(/nx_in, ny_in , 1, 1/), count=(/ngridx, ngridy, 1, ntime/) )
     case('PSFC')
        status = nf90_get_var(ncid, VarId, press_tmp(:, :, 0:0, :), &
             start=(/nx_in, ny_in , 1, 1/), count=(/ngridx, ngridy, 1, ntime/) )
     case('Q2')
        status = nf90_get_var(ncid, VarId, mixratio(:, :, 0:0, :), &
             start=(/nx_in, ny_in , 1, 1/), count=(/ngridx, ngridy, 1, ntime/) )
     case('XLAT')
        status = nf90_get_var(ncid, VarId, lat, &
             start=(/nx_in, ny_in, 1/), count=(/ngridx, ngridy, 1/))
     case('XLONG')
        status = nf90_get_var(ncid, VarId, lon, &
             start=(/nx_in, ny_in, 1/), count=(/ngridx, ngridy, 1/))
        ! *** Reading for WRF PROFILE variables:
     case('PHB')
        status = nf90_get_var(ncid, VarId, hgt_tmp, &
             start=(/nx_in, ny_in , 1, 1/), &
             count=(/ngridx, ngridy, nlyr+1, ntime/) )
        hgt_tmp = 1.0E-3*hgt_tmp/9.81  ! [km]

     case('P_HYD')
        status = nf90_get_var(ncid, VarId, press_tmp(:, :, 1:nlyr, :), &
             start=(/nx_in, ny_in , 1, 1/), count=(/ngridx, ngridy, nlyr, ntime/) )
        press_tmp = press_tmp*1.E-2  ! [hPa]
     case('T')
        status = nf90_get_var(ncid, VarId, pertur_T, &
             start=(/nx_in, ny_in , 1, 1/), count=(/ngridx, ngridy, nlyr, ntime/) )
     case('QVAPOR')
        status = nf90_get_var(ncid, VarId, mixratio(:, :, 1:nlyr, :), &
             start=(/nx_in, ny_in , 1, 1/), count=(/ngridx, ngridy, nlyr, ntime/) )
     case('QCLOUD')
        status = nf90_get_var(ncid, VarId, cloud_water_tmp, &
             start=(/nx_in, ny_in , 1, 1/), count=(/ngridx, ngridy, nlyr, ntime/) )
     case('QRAIN')
        status = nf90_get_var(ncid, VarId, rain_water_tmp, &
             start=(/nx_in, ny_in , 1, 1/), count=(/ngridx, ngridy, nlyr, ntime/) )
     case('QICE')
        status = nf90_get_var(ncid, VarId, cloud_ice_tmp, &
             start=(/nx_in, ny_in , 1, 1/), count=(/ngridx, ngridy, nlyr, ntime/) )
     case('QSNOW')
        status = nf90_get_var(ncid, VarId, snow_tmp, &
             start=(/nx_in, ny_in , 1, 1/), count=(/ngridx, ngridy, nlyr, ntime/) )
     case('QGRAUP')
        status = nf90_get_var(ncid, VarId, graupel_tmp, &
             start=(/nx_in, ny_in , 1, 1/), count=(/ngridx, ngridy, nlyr, ntime/) )
     case('V')
        status = nf90_get_var(ncid, VarId, V_Vel, &
             start=(/nx_in, ny_in+1 , 1, 1/), count=(/ngridx, ngridy, nlyr, ntime/) )
     case('U')
        status = nf90_get_var(ncid, VarId, U_Vel, &
             start=(/nx_in+1, ny_in , 1, 1/), count=(/ngridx, ngridy, nlyr, ntime/) )
     case('Times')
        status = nf90_get_var(ncid, VarId, TimeStamp)

     case default
        print*, 'WARNING: WRF variable ', trim(varname),' not recognized.'
     end select
  end do


  ! ---
  ! Variables not present in WRF:
  ! 1) Converting Perturbation Potential Temperature to Temperature:
  temp_tmp(:, :, 1:nlyr, :) = PERTHETA2T(pertur_T, &
       & press_tmp(:, :, 1:nlyr, :) )
  
  ! 2) Converting Vapor mixing ratio to Relative Humidity
  relhum_tmp = mixr_to_rh(mixratio, press_tmp, temp_tmp)

  ! 3) Converting Wind U and V components to Windspeed and Direction:
  call Wind_UV2speeddir(U_Vel, V_vel, windvel_tmp, winddir_tmp)

  ! 4) Converting TimeStamp to Vector Date variable:
  UnixTime = getUnixTime(TimeStamp)

  ! --
  ! 5) Converting mass-ratio [kg/kg] to mass-volume [g/m^3] values:
  mixr_tmp = mixratio(:, :, 1:nlyr, :)  ! [kg/kg]
  call Convert_Hydrometeors_units(Is_Specific_VAR=.false.)

  ! ****************************
  ! Retrieving Global Attribute:
  status = nf90_get_att(ncid,NF90_GLOBAL, 'DX', del_xy(1))
  status = nf90_get_att(ncid,NF90_GLOBAL, 'DY', del_xy(2))
  status = nf90_get_att(ncid,NF90_GLOBAL, 'TITLE', varname)
  write(origin_str,'(A)')  trim(varname)//'->'//trim(ncfile)
  
  ! ****************************
  ! Closing the NetCDF file
  status = nf90_close(ncid)
  call check_nc(status)
  
  deallocate(pertur_T, U_Vel, V_Vel, mixratio, TimeStamp)
  
  return
end subroutine read_wrf
! _____________________________________________________________________
! _____________________________________________________________________


! _____________________________________________________________________
! ---------------------------------------------------------------------
! Subroutine to read Wyoming Radiosonde database from the NetCDF-file generated by the
! code in repository github.com/pablosaa/WyoSondes
!
! The loaded variables are introduced to the RT3/RT4 radiative transfer code for
! calculation of Brightness Temperature Fields.
! ----------------------------------------------------------------------
subroutine read_wyosonde(ncflen,ncfile,mxgridx,mxgridy,mxlyr,mxtime,hgt_lev,&
     &press_lev,temp_lev,relhum_lev,mixr_lev, cloud_water_lev,&
     &rain_water_lev,cloud_ice_lev,snow_lev,graupel_lev, winddir_lev, windvel_lev,&
     &qidx, ngridx, ngridy, del_xy, nlyr, ntime, lat, lon, year, month, day, hour,&
     &origin_str)
  use netcdf
  implicit none

  integer, intent(in) :: ncflen
  character(len=ncflen), intent(in) :: ncfile

  ! -- Here come variables declaration alike the onses used by RT3

  integer, intent(in) ::   mxgridx, mxgridy, mxlyr, mxtime
  real(kind=8), intent(inout), dimension(mxgridx,mxgridy,0:mxlyr,mxtime) :: hgt_lev,temp_lev,press_lev,relhum_lev
  real(kind=8), intent(inout), dimension(mxgridx,mxgridy,mxlyr,mxtime) :: mixr_lev, cloud_water_lev,rain_water_lev,cloud_ice_lev
  real(kind=8), intent(inout), dimension(mxgridx,mxgridy,mxlyr,mxtime) :: snow_lev,graupel_lev, winddir_lev, windvel_lev
  integer(kind=4), intent(inout), dimension(mxgridx,mxgridy,mxtime) :: qidx
  integer, intent(out) :: ngridx, ngridy, nlyr, ntime
  real(kind=4), intent(out), dimension(mxgridx, mxgridy) :: lat, lon
  real(kind=4), intent(out), dimension(mxgridx, mxgridy, mxtime) :: year, month, day, hour
  real(kind=4), intent(out) :: del_xy(2)
  character(len=*), intent(out) :: origin_str
  ! -- end of variable declaration

  integer :: status
  integer :: ncid, ndims_in, nvars_in, ngatts_in, unlimdimid_in
  integer :: VarId
  character(len=10) :: dim_name(4)
  character(len=10) :: varname
  integer :: i, NN, dim_len(4)
  integer, allocatable, dimension(:) :: myVarIDs
  real, allocatable, dimension(:,:,:,:) :: Var4D
  real, allocatable, dimension(:,:,:) :: Var3D
  real, allocatable, dimension(:,:) :: Var2D
  real, allocatable, dimension(:) :: Var1D

  ! Initialazing the variables to a fixed value
  hgt_lev = 0.
  temp_lev =0.
  press_lev = 0.
  relhum_lev = 0.
  cloud_water_lev = 0.
  rain_water_lev = 0.
  cloud_ice_lev = 0.
  snow_lev = 0.
  graupel_lev = 0.
  lat = 0.
  lon = 0.
  year = 0.
  month = 0.
  day = 0.
  hour = 0.

  ! -----------------------------------------------------------------------------
  ! Open file and read directory
  print*,'netCDF input files is', ncflen, ' : ', ncfile
  
  status = nf90_open(ncfile,NF90_NOWRITE,ncid)
  if(status.ne.0) then
     print*, 'Cannot open netCDF-file: '//trim(ncfile)
     stop
  end if
  status = nf90_inquire(ncid,ndims_in,nvars_in,ngatts_in,unlimdimid_in)
  ! Get ID of unlimited dimension

  !! print*,ndims_in, nvars_in, ngatts_in, unlimdimid_in
  do i=1,ndims_in
     ! assigning dimension name and length:
     status = nf90_inquire_dimension(ncid,i,dim_name(i), dim_len(i))
     select case(trim(dim_name(i)))
     case('xn')
        ngridx = dim_len(i)
     case('yn')
        ngridy = dim_len(i)
     case('lev')
        nlyr = dim_len(i)
     case('time')
        ntime = dim_len(i)
     case default
        print*, 'netCDF file dimension '//trim(dim_name(i))//' unknown!'
        stop
     end select
  end do
  allocate(myVarIDs(nvars_in))
  status = nf90_inq_varids(ncid,nvars=NN,varids=myVarIDs)

  ! Reading variables:
  allocate(Var4D(dim_len(1),dim_len(2),dim_len(3),dim_len(4)))
  allocate(Var3D(dim_len(1),dim_len(2),dim_len(4)))
  allocate(Var2D(dim_len(1),dim_len(2)))
  allocate(Var1D(dim_len(4)))

  do i=1,nvars_in
     status = nf90_inquire_variable(ncid, myVarIDs(i), varname, ndims = NN)
     if(status /= nf90_NoErr) print*, 'ERROR: NetCDF variable name cannot be assigned'
     status = nf90_inq_varid(ncid, varname, VarId)


     if(status /= nf90_NoErr) print*, 'ERROR: NetCDF variable ID for ',varname,' cannot be retrieved'

     Var4D = -99.
     Var3D = -99.
     Var2D = -99.
     Var1D = -99.

     ! Loading the variables from NetCDF
     select case(NN)
     case(1)
        status = nf90_get_var(ncid, VarId, Var1D)
     case(2)
        status = nf90_get_var(ncid, VarId, Var2D)
     case(3)
        status = nf90_get_var(ncid, VarId, Var3D)
     case(4)
        status = nf90_get_var(ncid, VarId, Var4D)
     case default
        print*,trim(varname),'WARNING: neither 4D nor 3D nor 2D variable!!'
        cycle
     end select
     if(status /= nf90_NoErr) print*, 'error getting variable ', trim(varname)

     ! Assigning the data to RT3/4 variable names
     select case(trim(varname))
     case('HGT')
        hgt_lev(:dim_len(1),:dim_len(2),0,:dim_len(4)) = spread(Var2D,dim=3,ncopies=ntime)
        ! for 2D variables:
     case('T2')
        temp_lev(:dim_len(1),:dim_len(2),0,:dim_len(4)) = Var3D
     case('PSFC')
        press_lev(:dim_len(1),:dim_len(2),0,:dim_len(4)) = Var3D
     case('Q2')
        relhum_lev(:dim_len(1),:dim_len(2),0,:dim_len(4)) = Var3D
        ! for 3D variables:
     case('PHB')
        ! Passing values to RT3 variables:
        hgt_lev(:dim_len(1),:dim_len(2),1:dim_len(3),:dim_len(4)) = Var4D
     case('P')
        press_lev = 0
        ! Passing values to RT3 variables:
        press_lev(:dim_len(1),:dim_len(2),1:dim_len(3),:dim_len(4)) = Var4D
     case('T')
        temp_lev = 0
        ! Passing values to RT3 variables:
        temp_lev(:dim_len(1),:dim_len(2),1:dim_len(3),:dim_len(4)) = Var4D
     case('RH')
        relhum_lev = 0
        ! Passing values to RT3 variables:
        relhum_lev(:dim_len(1),:dim_len(2),1:dim_len(3),:dim_len(4)) = Var4D
     case('QVAPOR')
        mixr_lev = 0
        ! Passing values to RT3 variables:
        mixr_lev(:dim_len(1),:dim_len(2),1:dim_len(3),:dim_len(4)) = Var4D
     case('QCLOUD')
        cloud_water_lev = 0
        ! Passing values to RT3 variables:
        cloud_water_lev(:dim_len(1),:dim_len(2),1:dim_len(3),:dim_len(4)) = Var4D
     case('QRAIN')
        rain_water_lev = 0
        ! Passing values to RT3 variables:
        rain_water_lev(:dim_len(1),:dim_len(2),1:dim_len(3),:dim_len(4)) = Var4D
     case('QICE')
        cloud_ice_lev = 0
        ! Passing values to RT3 variables:
        cloud_ice_lev(:dim_len(1),:dim_len(2),1:dim_len(3),:dim_len(4)) = Var4D
     case('QSNOW')
        snow_lev = 0
        ! Passing values to RT3 variables:
        snow_lev(:dim_len(1),:dim_len(2),1:dim_len(3),:dim_len(4)) = Var4D
     case('QGRAUP')
        graupel_lev = 0
        ! Passing values to RT3 variables:
        graupel_lev(:dim_len(1),:dim_len(2),1:dim_len(3),:dim_len(4)) = Var4D
     case('WINDDIR')
        winddir_lev = 0
        ! Passing values to RT3 variables:
        winddir_lev(:dim_len(1),:dim_len(2),1:dim_len(3),:dim_len(4)) = Var4D
     case('WINDVEL')
        windvel_lev = 0
        ! Passing values to RT3 variables:
        windvel_lev(:dim_len(1),:dim_len(2),1:dim_len(3),:dim_len(4)) = Var4D
     case('QIDX')
        ! Passing values to quality index:
        qidx(:dim_len(1),:dim_len(2),:dim_len(4)) = Var3D
     case('LAT')
        lat(:dim_len(1),:dim_len(2)) = Var2D
     case('LON')
        lon(:dim_len(1),:dim_len(2)) = Var2D
     case('year')
        year(:dim_len(1),:dim_len(2),:dim_len(4)) = Var3D
     case('month')
        month(:dim_len(1),:dim_len(2),:dim_len(4)) = Var3D
     case('day')
        day(:dim_len(1),:dim_len(2),:dim_len(4)) = Var3D
     case('hour')
        hour(:dim_len(1),:dim_len(2),:dim_len(4)) = Var3D
     case default
        print*, 'variable ',varname,' not assigned yet'
     end select
  end do

  ! ****************************
  ! Retrieving Global Attribute:
  status = nf90_get_att(ncid,NF90_GLOBAL,'grid_x',del_xy(1))
  status = nf90_get_att(ncid,NF90_GLOBAL,'grid_y',del_xy(2))
  status = nf90_get_att(ncid,NF90_GLOBAL,'origin',origin_str)
  
  ! ****************************
  ! Closing the NetCDF file
  status = nf90_close(ncid)
  
  deallocate(myVarIDs, Var4D, Var3D, Var2D, Var1D)

  return
end subroutine read_wyosonde
! _____________________________________________________________________


! =====================================================================
! ---------------------------------------------------------------------
! SUBROUTINE to create the RT3/4 output as NetCDF files 
!
! This subroutine only creates the NetCDF file with its corresponding
! variable names, dimentision and attributes.
! When the main code needs to write data, then the subroutine named
! 'storencdf' needs to be called in order to pass the FORTRAN variables
! to their corresponding NetCDF variables.
!
! ---------------------------------------------------------------------
!subroutine createncdf(ncflen, ncfile,NUMMU,NFREQ,NSTOKES,NLYR,XN,YN,&
!     &LAYERS,freq_str,input_file,micro_phys,SELV,SLAT,SLON,origin_str)
subroutine createncdf(ncflen, ncfile,NUMMU,NSTOKES,&
     input_file,micro_phys,origin_str)
  use netcdf
  use nctoys, only : check_nc
  use variables, only : nelv, elevations, n_freq, nlyr, ngridx, ngridy, hgt_tmp,&
       FREQ, lat, lon, PBLH, qidx, temp_tmp, press_tmp, relhum_tmp, mixr_tmp, &
       cloud_water_tmp, rain_water_tmp, cloud_ice_tmp, snow_tmp, graupel_tmp, &
       winddir_tmp, windvel_tmp, ntime
  implicit none

  integer, intent(in) :: ncflen
  character(len=ncflen), intent(in) :: ncfile
  character(len=*) input_file, micro_phys
  integer, intent(in) :: NUMMU, NSTOKES
  character(len=*), intent(in) :: origin_str
  
  ! internal variables
  integer :: I, status, ncid
  integer :: nDims, mu_id, stok_id, freq_id, lyr_id, time_id, xn_id, yn_id
  integer :: var_mu_id, var_stok_id, var_freq_id, var_lyr_id, var_time_id, var_elvid, var_latid, var_lonid, var_gph_id
  integer :: var_xid, var_yid, var_tbup1_id, var_tbdn1_id, var_tbup0_id, var_tbdn0_id
  integer :: var_te2_id, var_rh2_id, var_pr2_id, var_blh_id
  integer :: var_te_id, var_pr_id, var_rh_id, var_qv_id, var_qc_id
  integer :: var_qr_id, var_qs_id, var_qg_id, var_qi_id, var_wd_id, var_ws_id
  integer :: var_kextqc_id, var_kextqr_id, var_kextqs_id
  integer :: var_kextqg_id, var_kextqi_id, var_kextatm_id
  integer :: var_salbtot_id, var_backsct_id, var_gcoeff_id, var_qidx_id
  integer, dimension(NSTOKES) :: stokes_var
  integer :: NANGLES
  character(len=8) :: ini_date
  character(len=10) :: ini_time

  NANGLES = NUMMU + nelv

  call date_and_time(DATE=ini_date, TIME=ini_time)
  status = nf90_create(trim(ncfile), ior(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid)
  if(status /= nf90_NOERR) stop 'Output NetCDF not possible to create: '  !//nf90_strerror(status) 

  ! Defining dimensions
  status = nf90_def_dim(ncid, "theta_z", NANGLES, mu_id) ! NUMMU
  status = nf90_def_dim(ncid, "freq", N_FREQ, freq_id)
  status = nf90_def_dim(ncid, "stokes", NSTOKES, stok_id)
  status = nf90_def_dim(ncid, "layer", NLYR, lyr_id)
  status = nf90_def_dim(ncid, "time", NF90_UNLIMITED, time_id)
  status = nf90_def_dim(ncid, "lat", ngridx , xn_id)
  status = nf90_def_dim(ncid, "lon", ngridy , yn_id)

  ! Define of variables
  status = nf90_def_var(ncid, "theta_z", NF90_REAL4, (/ mu_id /), var_mu_id)
  status = nf90_def_var(ncid, "freq", NF90_REAL4, (/ freq_id /), var_freq_id)
  status = nf90_def_var(ncid, "stokes", NF90_INT, (/ stok_id /), var_stok_id)
  status = nf90_def_var(ncid, "layer", NF90_REAL4, (/lyr_id/), var_lyr_id)
  status = nf90_def_var(ncid, "GPH", NF90_REAL, (/xn_id, yn_id, lyr_id, time_id/), var_gph_id)
  status = nf90_def_var(ncid, "time", NF90_REAL, (/ time_id /), var_time_id)
  status = nf90_def_var(ncid, "TB_UP_TOA", NF90_REAL4, (/ mu_id, freq_id, stok_id, xn_id, yn_id, time_id /), var_tbup1_id)
  status = nf90_def_var(ncid, "TB_UP_GRD", NF90_REAL4, (/ mu_id, freq_id, stok_id, xn_id, yn_id, time_id /), var_tbup0_id)
  status = nf90_def_var(ncid, "TB_DN_GRD", NF90_REAL4, (/ mu_id, freq_id, stok_id, xn_id, yn_id, time_id /), var_tbdn0_id)
  status = nf90_def_var(ncid, "TB_DN_TOA", NF90_REAL4, (/ mu_id, freq_id, stok_id, xn_id, yn_id, time_id /), var_tbdn1_id)

  !
  ! Define atmospheric state variables
  ! 1. Station level variables 
  status = nf90_def_var(ncid, "T2m", NF90_REAL4, (/xn_id, yn_id, time_id/), var_te2_id)
  status = nf90_def_var_deflate(ncid, var_te2_id, 0, 1, 9)
  status = nf90_def_var(ncid, "RH2m", NF90_REAL4, (/xn_id, yn_id, time_id/), var_rh2_id)
  status = nf90_def_var_deflate(ncid, var_rh2_id, 0, 1, 9)
  status = nf90_def_var(ncid, "P2m", NF90_REAL4, (/xn_id, yn_id, time_id/), var_pr2_id)
  status = nf90_def_var_deflate(ncid, var_pr2_id, 0, 1, 9)
  status = nf90_def_var(ncid, "PBLH", NF90_REAL4, (/xn_id, yn_id, time_id/), var_blh_id)
  status = nf90_def_var_deflate(ncid, var_blh_id, 0, 1, 9)
  
  ! 2. Atmospheric Profile variables 
  status = nf90_def_var(ncid, "temp", NF90_REAL4, (/xn_id, yn_id, lyr_id, time_id/), var_te_id)
  status = nf90_def_var_deflate(ncid, var_te_id, 0, 1, 9)
  status = nf90_def_var(ncid, "press", NF90_REAL4, (/xn_id, yn_id, lyr_id, time_id/), var_pr_id)
  status = nf90_def_var_deflate(ncid, var_pr_id, 0, 1, 9)
  status = nf90_def_var(ncid, "rh", NF90_REAL4, (/xn_id, yn_id, lyr_id, time_id/), var_rh_id)
  status = nf90_def_var_deflate(ncid, var_rh_id, 0, 1, 9)
  status = nf90_def_var(ncid, "qv", NF90_REAL4, (/xn_id, yn_id, lyr_id, time_id/), var_qv_id)
  status = nf90_def_var_deflate(ncid, var_qv_id, 0, 1, 9)
  status = nf90_def_var(ncid, "qc", NF90_REAL4, (/xn_id, yn_id, lyr_id, time_id/), var_qc_id)
  status = nf90_def_var_deflate(ncid, var_qc_id, 0, 1, 9)
  status = nf90_def_var(ncid, "qr", NF90_REAL4, (/xn_id, yn_id, lyr_id, time_id/), var_qr_id)
  status = nf90_def_var_deflate(ncid, var_qr_id, 0, 1, 9)
  status = nf90_def_var(ncid, "qs", NF90_REAL4, (/xn_id, yn_id, lyr_id, time_id/), var_qs_id)
  status = nf90_def_var_deflate(ncid, var_qs_id, 0, 1, 9)
  status = nf90_def_var(ncid, "qg", NF90_REAL4, (/xn_id, yn_id, lyr_id, time_id/), var_qg_id)
  status = nf90_def_var_deflate(ncid, var_qg_id, 0, 1, 9)
  status = nf90_def_var(ncid, "qi", NF90_REAL4, (/xn_id, yn_id, lyr_id, time_id/), var_qi_id)
  status = nf90_def_var_deflate(ncid, var_qi_id, 0, 1, 9)
  status = nf90_def_var(ncid, "wd", NF90_REAL4, (/xn_id, yn_id, lyr_id, time_id/), var_wd_id)
  status = nf90_def_var_deflate(ncid, var_wd_id, 0, 1, 9)
  status = nf90_def_var(ncid, "ws", NF90_REAL4, (/xn_id, yn_id, lyr_id, time_id/), var_ws_id)
  status = nf90_def_var_deflate(ncid, var_ws_id, 0, 1, 9)
  
  ! Definition of Micro-physics frequency dependent variables (Profiles)
  status = nf90_def_var(ncid, "kext_atm", NF90_REAL4, (/xn_id, yn_id, lyr_id, freq_id, time_id/), var_kextatm_id)
  status = nf90_def_var_deflate(ncid, var_kextatm_id, 0, 1, 9)
  status = nf90_def_var(ncid, "kext_qc", NF90_REAL4, (/xn_id, yn_id, lyr_id, freq_id, time_id/), var_kextqc_id)
  status = nf90_def_var_deflate(ncid, var_kextqc_id, 0, 1, 9)
  status = nf90_def_var(ncid, "kext_qr", NF90_REAL4, (/xn_id, yn_id, lyr_id, freq_id, time_id/), var_kextqr_id)
  status = nf90_def_var_deflate(ncid, var_kextqr_id, 0, 1, 9)
  status = nf90_def_var(ncid, "kext_qs", NF90_REAL4, (/xn_id, yn_id, lyr_id, freq_id, time_id/), var_kextqs_id)
  status = nf90_def_var_deflate(ncid, var_kextqs_id, 0, 1, 9)
  status = nf90_def_var(ncid, "kext_qg", NF90_REAL4, (/xn_id, yn_id, lyr_id, freq_id, time_id/), var_kextqg_id)
  status = nf90_def_var_deflate(ncid, var_kextqg_id, 0, 1, 9)
  status = nf90_def_var(ncid, "kext_qi", NF90_REAL4, (/xn_id, yn_id, lyr_id, freq_id, time_id/), var_kextqi_id)
  status = nf90_def_var_deflate(ncid, var_kextqi_id, 0, 1, 9)
  status = nf90_def_var(ncid, "alb_tot", NF90_REAL4, (/xn_id, yn_id, lyr_id, freq_id, time_id/), var_salbtot_id)
  status = nf90_def_var_deflate(ncid, var_salbtot_id, 0, 1, 9)
  status = nf90_def_var(ncid, "back_scatt", NF90_REAL4, (/xn_id, yn_id, lyr_id, freq_id, time_id/), var_backsct_id)
  status = nf90_def_var_deflate(ncid, var_backsct_id, 0, 1, 9)
  status = nf90_def_var(ncid, "g_coeff", NF90_REAL4, (/xn_id, yn_id, lyr_id, freq_id, time_id/), var_gcoeff_id)
  status = nf90_def_var_deflate(ncid, var_gcoeff_id, 0, 1, 9)

  ! profile quality index
  status = nf90_def_var(ncid, "QIDX", NF90_INT4, (/xn_id, yn_id, time_id/), var_qidx_id)
  status = nf90_def_var_deflate(ncid, var_qidx_id, 0, 1, 9)
  
  ! grid-based variables
  !status = nf90_def_var(ncid, "xn", NF90_INT, (/ xn_id /), var_xid)  ! time_id
  !status = nf90_def_var(ncid, "yn", NF90_INT, (/ yn_id /), var_yid)   ! time_id

  status = nf90_def_var(ncid, "elevation", NF90_REAL4, (/xn_id, yn_id/), var_elvid)
  status = nf90_def_var_deflate(ncid, var_elvid, 0, 1, 9)
  status = nf90_def_var(ncid, "lat", NF90_REAL4, (/xn_id, yn_id/), var_latid)
  status = nf90_def_var_deflate(ncid, var_latid, 0, 1, 9)
  status = nf90_def_var(ncid, "lon", NF90_REAL4, (/xn_id, yn_id/), var_lonid)
  status = nf90_def_var_deflate(ncid, var_lonid, 0, 1, 9)

  ! ***********************************************
  ! Adding Attributes for Variables
  status = nf90_put_att(ncid, var_mu_id, "short_name","theta_z")
  status = nf90_put_att(ncid, var_mu_id, "long_name","Zenithal angle")
  status = nf90_put_att(ncid, var_mu_id, "units","degree")
  status = nf90_put_att(ncid, var_mu_id, "N_obs_angles", nelv)
  status = nf90_put_att(ncid, var_mu_id, "Obs_angles_degree", elevations(:nelv))
  status = nf90_put_att(ncid, var_mu_id, "N_sim_angles",NUMMU)

  status = nf90_put_att(ncid, var_stok_id, "short_name","stk")
  status = nf90_put_att(ncid, var_stok_id, "long_name","Stokes_vector")
  status = nf90_put_att(ncid, var_stok_id, "units","1")

  status = nf90_put_att(ncid, var_lyr_id, "short_name","Layers")
  status = nf90_put_att(ncid, var_lyr_id, "long_name"," Layer height at center grid")
  status = nf90_put_att(ncid, var_lyr_id, "units","km")
  status = nf90_put_att(ncid, var_lyr_id, "_FillValue", -999.9)
  status = nf90_put_att(ncid, var_lyr_id, "Note", "For reference only")

  status = nf90_put_att(ncid, var_gph_id, "short_name","h.a.g.l")
  status = nf90_put_att(ncid, var_gph_id, "long_name","Geopotential height")
  status = nf90_put_att(ncid, var_gph_id, "units","km")
  status = nf90_put_att(ncid, var_gph_id, "_FillValue", -999.9)

  status = nf90_put_att(ncid, var_freq_id, "short_name","freq")
  status = nf90_put_att(ncid, var_freq_id, "long_name","Radiometric frequency")
  status = nf90_put_att(ncid, var_freq_id, "units","GHz")

  status = nf90_put_att(ncid, var_time_id, "short_name","time")
  status = nf90_put_att(ncid, var_time_id, "long_name","days since 1970.1.1 00:00:00")
  status = nf90_put_att(ncid, var_time_id, "units","day")
  status = nf90_put_att(ncid, var_time_id, "_FillValue", -999.9)

  status = nf90_put_att(ncid, var_tbup1_id, "short_name","TB_UP_TOA")
  status = nf90_put_att(ncid, var_tbup1_id, "long_name","TOA Brightness Temperature (Upwelling)")
  status = nf90_put_att(ncid, var_tbup1_id, "units","K")
  status = nf90_put_att(ncid, var_tbup1_id, "_FillValue", -999.9)

  status = nf90_put_att(ncid, var_tbdn1_id, "short_name","TB_DN_TOA")
  status = nf90_put_att(ncid, var_tbdn1_id, "long_name","TOA Brightness Temperature (Downwelling)")
  status = nf90_put_att(ncid, var_tbdn1_id, "units","K")
  status = nf90_put_att(ncid, var_tbdn1_id, "_FillValue", -999.9)

  status = nf90_put_att(ncid, var_tbup0_id, "short_name","TB_UP_GND")
  status = nf90_put_att(ncid, var_tbup0_id, "long_name","GROUND Brightness Temperature (Upwelling)")
  status = nf90_put_att(ncid, var_tbup0_id, "units","K")
  status = nf90_put_att(ncid, var_tbup0_id, "_FillValue", -999.9)

  status = nf90_put_att(ncid, var_tbdn0_id, "short_name","TB_DN_GND")
  status = nf90_put_att(ncid, var_tbdn0_id, "long_name","GROUND Brightness Temperature (Downwelling)")
  status = nf90_put_att(ncid, var_tbdn0_id, "units","K")
  status = nf90_put_att(ncid, var_tbdn0_id, "_FillValue", -999.9)

  ! Attributes for Station level variables
  status = nf90_put_att(ncid, var_te2_id, "short_name","T2m")
  status = nf90_put_att(ncid, var_te2_id, "long_name","2m Temperature")
  status = nf90_put_att(ncid, var_te2_id, "units","K")
  status = nf90_put_att(ncid, var_te2_id, "_FillValue",-999.9)

  status = nf90_put_att(ncid, var_rh2_id, "short_name","RH_s")
  status = nf90_put_att(ncid, var_rh2_id, "long_name","2m Relative humidity")
  status = nf90_put_att(ncid, var_rh2_id, "units","%")
  status = nf90_put_att(ncid, var_rh2_id, "_FillValue",-999.9)

  status = nf90_put_att(ncid, var_pr2_id, "short_name","P_s")
  status = nf90_put_att(ncid, var_pr2_id, "long_name","2m air pressure")
  status = nf90_put_att(ncid, var_pr2_id, "units","hPa")
  status = nf90_put_att(ncid, var_pr2_id, "_FillValue",-999.9)

  status = nf90_put_att(ncid, var_blh_id, "short_name","PBLH")
  status = nf90_put_att(ncid, var_blh_id, "long_name","Height of the PBL")
  status = nf90_put_att(ncid, var_blh_id, "units","m")
  status = nf90_put_att(ncid, var_blh_id, "_FillValue",-999.9)

  
  ! Attributes for Profile variables
  status = nf90_put_att(ncid, var_te_id, "short_name","temp")
  status = nf90_put_att(ncid, var_te_id, "long_name","temperature")
  status = nf90_put_att(ncid, var_te_id, "units","K")
  status = nf90_put_att(ncid, var_te_id, "_FillValue",-999.9)

  status = nf90_put_att(ncid, var_pr_id, "short_name","press")
  status = nf90_put_att(ncid, var_pr_id, "long_name","pressure")
  status = nf90_put_att(ncid, var_pr_id, "units","hPa")
  status = nf90_put_att(ncid, var_pr_id, "_FillValue",-999.9)

  status = nf90_put_att(ncid, var_rh_id, "short_name","rh")
  status = nf90_put_att(ncid, var_rh_id, "long_name","relative humidity")
  status = nf90_put_att(ncid, var_rh_id, "units","%")
  status = nf90_put_att(ncid, var_rh_id, "_FillValue",-999.9)
  
  status = nf90_put_att(ncid, var_qv_id, "short_name","qv")
  status = nf90_put_att(ncid, var_qv_id, "long_name","vapour mixing ratio")
  status = nf90_put_att(ncid, var_qv_id, "units","kg/kg")
  status = nf90_put_att(ncid, var_qv_id, "_FillValue",-999.9)

  status = nf90_put_att(ncid, var_qc_id, "short_name","qc")
  status = nf90_put_att(ncid, var_qc_id, "long_name","cloud water content")
  status = nf90_put_att(ncid, var_qc_id, "units","g m-3")
  status = nf90_put_att(ncid, var_qc_id, "_FillValue",-999.9)

  status = nf90_put_att(ncid, var_qr_id, "short_name","qr")
  status = nf90_put_att(ncid, var_qr_id, "long_name","rain water content")
  status = nf90_put_att(ncid, var_qr_id, "units","g m-3")
  status = nf90_put_att(ncid, var_qr_id, "_FillValue",-999.9)
  
  status = nf90_put_att(ncid, var_qs_id, "short_name","qs")
  status = nf90_put_att(ncid, var_qs_id, "long_name","snow water content")
  status = nf90_put_att(ncid, var_qs_id, "units","g m-3")
  status = nf90_put_att(ncid, var_qs_id, "_FillValue",-999.9)

  status = nf90_put_att(ncid, var_qg_id, "short_name","qg")
  status = nf90_put_att(ncid, var_qg_id, "long_name","graupel water content")
  status = nf90_put_att(ncid, var_qg_id, "units","g m-3")
  status = nf90_put_att(ncid, var_qg_id, "_FillValue",-999.9)

  status = nf90_put_att(ncid, var_qi_id, "short_name","qi")
  status = nf90_put_att(ncid, var_qi_id, "long_name","ice water content")
  status = nf90_put_att(ncid, var_qi_id, "units","g m-3")
  status = nf90_put_att(ncid, var_qi_id, "_FillValue",-999.9)

  status = nf90_put_att(ncid, var_wd_id, "short_name","wd")
  status = nf90_put_att(ncid, var_wd_id, "long_name","wind direction")
  status = nf90_put_att(ncid, var_wd_id, "units","deg")
  status = nf90_put_att(ncid, var_wd_id, "_FillValue",-999.9)

  status = nf90_put_att(ncid, var_ws_id, "short_name","ws")
  status = nf90_put_att(ncid, var_ws_id, "long_name","wind speed")
  status = nf90_put_att(ncid, var_ws_id, "units","knot")
  status = nf90_put_att(ncid, var_ws_id, "_FillValue",-999.9)
  
  status = nf90_put_att(ncid, var_kextatm_id, "short_name","kext_atm")
  status = nf90_put_att(ncid, var_kextatm_id, "long_name","atmospheric extinction coefficient")
  status = nf90_put_att(ncid, var_kextatm_id, "units","km-1")
  status = nf90_put_att(ncid, var_kextatm_id, "_FillValue",-999.9)

  status = nf90_put_att(ncid, var_kextqc_id, "short_name","kext_cloud")
  status = nf90_put_att(ncid, var_kextqc_id, "long_name","cloud extinction coefficient")
  status = nf90_put_att(ncid, var_kextqc_id, "units","km-1")
  status = nf90_put_att(ncid, var_kextqc_id, "_FillValue",-999.9)

  status = nf90_put_att(ncid, var_kextqr_id, "short_name","kext_rain")
  status = nf90_put_att(ncid, var_kextqr_id, "long_name","rain extinction coefficient")
  status = nf90_put_att(ncid, var_kextqr_id, "units","km-1")
  status = nf90_put_att(ncid, var_kextqr_id, "_FillValue",-999.9)
  
  status = nf90_put_att(ncid, var_kextqs_id, "short_name","kext_snow")
  status = nf90_put_att(ncid, var_kextqs_id, "long_name","snow extinction coefficient")
  status = nf90_put_att(ncid, var_kextqs_id, "units","km-1")
  status = nf90_put_att(ncid, var_kextqs_id, "_FillValue",-999.9)

  status = nf90_put_att(ncid, var_kextqg_id, "short_name","kext_graupel")
  status = nf90_put_att(ncid, var_kextqg_id, "long_name","graupel extinction coefficient")
  status = nf90_put_att(ncid, var_kextqg_id, "units","km-1")
  status = nf90_put_att(ncid, var_kextqg_id, "_FillValue",-999.9)

  status = nf90_put_att(ncid, var_kextqi_id, "short_name","kext_ice")
  status = nf90_put_att(ncid, var_kextqi_id, "long_name","ice extinction coefficient")
  status = nf90_put_att(ncid, var_kextqi_id, "units","km-1")
  status = nf90_put_att(ncid, var_kextqi_id, "_FillValue",-999.9)
  
  status = nf90_put_att(ncid, var_salbtot_id, "short_name","alb_tot")
  status = nf90_put_att(ncid, var_salbtot_id, "long_name","total surface albedo")
  status = nf90_put_att(ncid, var_salbtot_id, "units","-")
  status = nf90_put_att(ncid, var_salbtot_id, "_FillValue",-999.9)

  status = nf90_put_att(ncid, var_backsct_id, "short_name","backscatt")
  status = nf90_put_att(ncid, var_backsct_id, "long_name","backscattering coefficient")
  status = nf90_put_att(ncid, var_backsct_id, "units","km-1")
  status = nf90_put_att(ncid, var_backsct_id, "_FillValue",-999.9)

  status = nf90_put_att(ncid, var_gcoeff_id, "short_name","g_coeff")
  status = nf90_put_att(ncid, var_gcoeff_id, "long_name","asymetry factor")
  status = nf90_put_att(ncid, var_gcoeff_id, "units", "-")
  status = nf90_put_att(ncid, var_gcoeff_id, "_FillValue",-999.9)

  ! ***** Quality index
  status = nf90_put_att(ncid, var_qidx_id, "short_name","qidx")
  status = nf90_put_att(ncid, var_qidx_id, "long_name","Quality index")
  status = nf90_put_att(ncid, var_qidx_id, "units", "-")
  status = nf90_put_att(ncid, var_qidx_id, "_FillValue",-99)
  status = nf90_put_att(ncid, var_qidx_id, "note", "15=full quality")

  ! ***** Grid variables 
  !status = nf90_put_att(ncid, var_xid, "short_name", "Xn")
  !status = nf90_put_att(ncid, var_xid, "long_name", "X_grid")
  !status = nf90_put_att(ncid, var_yid, "short_name", "Yn")
  !status = nf90_put_att(ncid, var_yid, "long_name", "Y_grid")

  status = nf90_put_att(ncid, var_elvid, "short_name", "ELV")
  status = nf90_put_att(ncid, var_elvid, "long_name", "Elevation a.s.l.")
  status = nf90_put_att(ncid, var_elvid, "units", "km")

  status = nf90_put_att(ncid, var_latid, "short_name", "LAT")
  status = nf90_put_att(ncid, var_latid, "long_name", "Latitude")
  status = nf90_put_att(ncid, var_latid, "units", "degree_north")
  
  status = nf90_put_att(ncid, var_lonid, "short_name", "LON")
  status = nf90_put_att(ncid, var_lonid, "long_name", "Longitude")
  status = nf90_put_att(ncid, var_lonid, "units", "degree_east")

  status = nf90_put_att(ncid,NF90_GLOBAL,"Input_data", input_file)
  status = nf90_put_att(ncid,NF90_GLOBAL,"Origin_Information", trim(origin_str))
  status = nf90_put_att(ncid,NF90_GLOBAL,"Hydrometeor_Microphysics", micro_phys)
  status = nf90_put_att(ncid,NF90_GLOBAL,"Contact","Pablo.Saavedra@uib.no")
  status = nf90_put_att(ncid,NF90_GLOBAL,"Institution","Geophysical Institute, University of Bergen")
  status = nf90_put_att(ncid,NF90_GLOBAL,"License","CC-BY-SA")

  ! writting Initial date as global variable:
  status = nf90_put_att(ncid,NF90_GLOBAL,"Start_Time", &
       ini_date//' '//ini_time )
  ! ********** End definitions *****************************
  status = nf90_enddef(ncid)
  ! ********************************************************
  ! Writting variables that won't change during calculation:
  status = nf90_put_var(ncid, var_qidx_id, qidx)
  status = nf90_put_var(ncid, var_blh_id, PBLH)
  

  ! Putting variables independent of time:
  stokes_var = (/(I,I=1,NSTOKES)/)
  !TIMELINE = -999.9
  status = nf90_put_var(ncid, var_lyr_id, &
       sum(hgt_tmp(int(ngridx/2),int(ngridy/2),1:nlyr,:), dim=2)/ntime )
  status = nf90_put_var(ncid, var_freq_id, FREQ)
  status = nf90_put_var(ncid, var_stok_id, stokes_var)
  !status = nf90_put_var(ncid, var_time_id, TIMELINE)
  status = nf90_put_var(ncid, var_elvid, real(hgt_tmp(:,:,0:0,1:1), 4) )
  status = nf90_put_var(ncid, var_latid, lat)
  status = nf90_put_var(ncid, var_lonid, lon)

  ! ********** Writting Station level Variables ****************************
  ! Geopotential Height
  status = nf90_put_var(ncid, var_gph_id, hgt_tmp(:,:,1:nlyr,:) )
  call check_nc(status, 'Geopotential cannot be written!', .true.)
  
  ! Writting the 2m Temperature
  status = nf90_put_var(ncid, var_te2_id, temp_tmp(:,:,0:0,:) )
  call check_nc(status, 'T2m variable cannot be written!', .true.)

  ! Writting the 2m RH
  status = nf90_put_var(ncid, var_rh2_id, relhum_tmp(:,:,0:0,:) )
  call check_nc(status, 'RH2m variable cannot be written!', .true.)
  
  ! Writting the 2m Air Pressure
  status = nf90_put_var(ncid, var_pr2_id, press_tmp(:,:,0:0,:) )
  call check_nc(status, 'P2m variable cannot be written!', .true.)
  
  ! ********** Writting Profile Variables ****************************
  ! Writting the Temperature
  status = nf90_put_var(ncid, var_te_id, temp_tmp(:,:,1:,:) )
  if(status /= nf90_NoErr) stop 'Temperature profile cannot be written!'
  
  ! Writting the Pressure
  status = nf90_put_var(ncid, var_pr_id, press_tmp(:,:,1:,:) )
  if(status /= nf90_NoErr) stop 'Pressure profile cannot be written!'

  ! Writting the Relative Humidity
  status = nf90_put_var(ncid, var_rh_id, relhum_tmp(:,:,1:,:) )
  if(status /= nf90_NoErr) stop 'RH profile cannot be written!'

  ! Writting the QV variable
  status = nf90_put_var(ncid, var_qv_id, mixr_tmp )  
  if(status /= nf90_NoErr) stop 'QV cannot be written!'
  
  ! Writting the QC variable
  status = nf90_put_var(ncid, var_qc_id, cloud_water_tmp )
  if(status /= nf90_NoErr) stop 'QC cannot be written!'

  ! Writting the QR variable
  status = nf90_put_var(ncid, var_qr_id, rain_water_tmp )
  if(status /= nf90_NoErr) stop 'QR cannot be written!'

  ! Writting the QI variable
  status = nf90_put_var(ncid, var_qi_id, cloud_ice_tmp )
  if(status /= nf90_NoErr) stop 'QI cannot be written!'

  ! Writting the QS variable
  status = nf90_put_var(ncid, var_qs_id, snow_tmp )
  if(status /= nf90_NoErr) stop 'QS cannot be written!'

  ! Writting the QG variable
  status = nf90_put_var(ncid, var_qg_id, graupel_tmp )
  if(status /= nf90_NoErr) stop 'QG cannot be written!'

  ! Writting the Wind Direction variable
  status = nf90_put_var(ncid, var_wd_id, winddir_tmp )
  if(status /= nf90_NoErr) stop 'WD cannot be written!'
  ! Writting the Wind Speed variable
  status = nf90_put_var(ncid, var_ws_id, windvel_tmp )
  if(status /= nf90_NoErr) stop 'WS cannot be written!'
  ! ---/
  
  status = nf90_close(ncid)
  call check_nc(status, 'Error closing after creation NetCDF', .true.)

end subroutine createncdf
! _____________________________________________________________________


! _____________________________________________________________________
! ---------------------------------------------------------------------
! SUBROUTINE to store the RT3/4 data in the created NetCDF output file
!
! ---------------------------------------------------------------------
subroutine storencdf(OUT_FILE,MU_VALUES,NUMMU,HEIGHT,NOUTLEVELS,OUTVAR,NSTOKES,time_len)
  use netcdf
  !use mpi
  use variables, only : nx_in, nx_fin, ny_in, ny_fin, timeidx, nx, ny, UnixTime
  use nctoys
  use meteo_tools, only : interp1_1D
  
  character(len=*), intent(in) :: OUT_FILE
  real(kind=8), intent(in) :: MU_VALUES(NUMMU), HEIGHT(NOUTLEVELS), OUTVAR(NUMMU,NSTOKES,2,NOUTLEVELS)
  integer, intent(in) :: NUMMU, NOUTLEVELS, NSTOKES, time_len

  ! internal variables
  integer :: i, j, k
  integer :: status, ncid, VarId, idx, idf
  integer :: nDims, unlimdimid
  character(len=len(OUT_FILE)+3) :: ncfile
  character(len=40) :: dim_name
  integer :: x_grid, y_grid, i_freq, NTIME
  real(kind=8), allocatable, dimension(:) :: ZENITH_THTA  ! (NUMMU)
  real(kind=8), allocatable, dimension(:,:,:,:) :: TB_THTA
  real(kind=8), allocatable, dimension(:) :: elevations, elvmu
  real(kind=8), parameter :: PI = 4.0*atan(1.0)
  integer :: nelv, NANG
  integer :: err
  
  ! Extracting information from the OUT_FILE character string:
  ! * The OUT_FILE has the form like:
  ! ../output/TB/RT3TB13090112Exp7.6MaxGa0.2Exp4.0MaxGaExp8.0x001y001f27.20
  !
  ! which needs to be transformed to a NetCDF file like:
  ! ../output/TB/RT3TB_Exp7.6MaxGa0.2Exp4.0MaxGaExp8.0_f27.2.nc

  !call MPI_Comm_size(MPI_COMM_WORLD, numprocs, err)
  !call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)

  ! * grid indices:
  idx = scan(OUT_FILE,'x',back=.true.)+1
  if(idx.eq.1) stop 'no x000_ found in string passed'

  read(OUT_FILE(idx:),'(I03XI03)') x_grid, y_grid
  if(x_grid.NE.nx.OR.y_grid.NE.ny) stop 'ERROR passing x_grid or y_grid in storecdf'
  
  !! idx = scan(OUT_FILE,'=',back=.true.)+1
  !! if(idx.eq.1) stop 'no separator = found in string passed'
  !! idf = scan(OUT_FILE,'x',back=.true.)
  !! read(OUT_FILE(idx:idf-1),'(4I02)') date(1,1:4) !, micro_phys  ! '(5X4I02A)'
  !! date(1,1) = date(1,1)+2000
  !! TIMELINE = F2UnixTime(date, ntime)

  ! constructing netCDF file to write data:
  idf = scan(OUT_FILE,'=',back=.true.)-1
  ncfile = OUT_FILE(:idf)

  status = NF90_OPEN(ncfile,MODE=IOR(NF90_WRITE,NF90_NETCDF4),NCID=ncid)
  !!!status = NF90_OPEN_PAR(ncfile, IOR(IOR(NF90_WRITE, NF90_NETCDF4), &
  !!!     NF90_MPIIO), MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
  call check_nc(status, 'STORENCDF() Opening the NetCDF', .true.)

  ! Getting NetCDF file dimensions and lengths:
  ! For time:
  !!!status = nf90_inquire(ncid, nDimensions = nDims,unlimitedDimID = unlimdimid)
  status = nf90_inq_varid(ncid,"time",unlimdimid)
  status = nf90_inquire_dimension(ncid,unlimdimid,dim_name, NTIME)

  ! For frequency
  ! * Frequency from OUT_FILE:
  idx = scan(OUT_FILE,'f',back=.true.)+1
  read(OUT_FILE(idx:),'(I6)') i_freq

  if(i_freq.LT.1) stop 'ERROR finding frequency index in STORENCDF'
  
  ! Writting variable values ZENITH_THTA into NetCDF file:
  status = nf90_inq_varid(ncid, "theta_z", VarId)
  call check_nc(status, 'cos(MU) ID cannot be assigned!')
  status = nf90_inquire_dimension(ncid,VarId,dim_name,NANG)
  status = nf90_get_att(ncid,VarID,"N_obs_angles",nelv)

  allocate(ZENITH_THTA(NANG))
  allocate(TB_THTA(NANG,NSTOKES,2,NOUTLEVELS))
  allocate(elevations(nelv))
  allocate(elvmu(nelv))

  status = nf90_get_att(ncid, VarId, "Obs_angles_degree", elevations)
  
  ! converting input MWR elevation angles into zenithal angles cos(pi/2-mu)
  elvmu = cos((90.0-elevations)*PI/180.)

  ! Checking whether additional observation angles are needed:
  if(nelv.GT.0) then
     ! Interpolate the values for the additional angles:
     CALL interp1_1D(MU_VALUES, OUTVAR, NUMMU, elvmu, nelv,&
          & ZENITH_THTA, NANG, TB_THTA, NSTOKES,2,NOUTLEVELS)
  else
     ! No additional angles, only simulations:
     ZENITH_THTA(:NANG) = MU_VALUES(:NUMMU)
  end if
  
  if(time_len.EQ.1.AND.i_freq.EQ.1) then
     ! converting cos(mu) to zenithal angle:
     ZENITH_THTA = acos(ZENITH_THTA)*180.0/PI
     status = nf90_var_par_access(ncid, VarId, NF90_INDEPENDENT)
     status = nf90_put_var(ncid,VarId,ZENITH_THTA(:NANG))
     call check_nc(status, 'cos(MU) values cannot be stored!', .true.)

  end if
    
  ! writing time
  ! writting TB_UPwelling
  status = nf90_inq_varid(ncid, "time", VarId)
  call check_nc(status, 'Time cannot be assigned!', .true.)
  status = nf90_var_par_access(ncid, VarId, NF90_INDEPENDENT)
  status = nf90_put_var(ncid,VarId, UnixTime(time_len), start=(/time_len/))
  call check_nc(status, 'Time values cannot be stored!', .true.)

  ! writting TOA TB_UPwelling
  ! OUTVAR has the dimension of [mu,stokes,Down/Up-welling, Obs_level]
  status = nf90_inq_varid(ncid, "TB_UP_TOA", VarId)
  call check_nc(status, 'TB_UP TOA cannot be assigned!', .true.)

  status = nf90_var_par_access(ncid, VarId, NF90_INDEPENDENT)
  call check_nc(status, 'TB_UP TOA cannot set collective mode')

  status = nf90_put_var(ncid, VarId, TB_THTA(:NANG,:,1:1,1:1), &
       start=(/1,i_freq,1,x_grid,y_grid,time_len/), &
       count=(/NANG,1,2,1,1,1/) )
  call check_nc(status, 'TB_UP TOA values cannot be stored!', .true.)
  
  ! writting TOA TB_DOWNwelling
  status = nf90_inq_varid(ncid, "TB_DN_TOA", VarId)
  call check_nc(status, 'TB_DN TOA cannot be assigned!', .true.)

  status = nf90_var_par_access(ncid, VarId, NF90_INDEPENDENT)
  call check_nc(status, 'TB_DN TOA cannot set collective mode')

  status = nf90_put_var(ncid, VarId, TB_THTA(:NANG,:,2:2,1:1), &
       start=(/1,i_freq,1,x_grid,y_grid,time_len/), &
       count=(/NANG,1,2,1,1,1/))
  call check_nc(status, 'TB_DN TOA values cannot be stored!', .true.)
  
  ! writting GROUND TB_UPwelling
  ! OUTVAR has the dimension [mu,stokes,Down/Up-welling,Obs_level]) 
  status = nf90_inq_varid(ncid, "TB_UP_GRD", VarId)
  call check_nc(status, 'TB_UP GRD cannot be assigned!', .true.)

  status = nf90_var_par_access(ncid, VarId, NF90_INDEPENDENT)
  call check_nc(status, 'TB_UP GRD cannot set collective mode')

  status = nf90_put_var(ncid, VarId, TB_THTA(:NANG,:,1:1,2:2), &
       start=(/1,i_freq,1,x_grid,y_grid,time_len/), &
       count=(/NANG,1,2,1,1,1/))
  call check_nc(status, 'TB_UP GRD values cannot be stored!', .true.)
  
  ! writting GROUND TB_DOWNwelling
  status = nf90_inq_varid(ncid, "TB_DN_GRD", VarId)
  call check_nc(status, 'TB_DN GRD cannot be assigned!', .true.)

  status = nf90_var_par_access(ncid, VarId, NF90_INDEPENDENT)
  call check_nc(status, 'TB_DN GRD cannot set collective mode')

  status = nf90_put_var(ncid, VarId, TB_THTA(:NANG,:, 2:2, 2:2), &
       start=(/1, i_freq, 1, x_grid, y_grid, time_len/), &
       count=(/NANG, 1, 2, 1, 1, 1/))
  call check_nc(status, 'TB_DN GRD values cannot be stored!', .true.)
    
  status = NF90_CLOSE(ncid)
  call check_nc(status, 'Closing (storencdf) not possible!', .true.)


  if(allocated(ZENITH_THTA)) deallocate(ZENITH_THTA)
  if(allocated(TB_THTA)) deallocate(TB_THTA)
  if(allocated(elevations)) deallocate(elevations)
  if(allocated(elvmu)) deallocate(elvmu)

  
  return
  
end subroutine storencdf
! ----/

! _____________________________________________________________________
! ---------------------------------------------------------------------
! SUBROUTINE Microphysics variable storege for the RT3/4 NetCDF output file
!
subroutine MP_storencdf(OUT_FILE, time_len, i_freq)
  use netcdf
  use nctoys, only : check_nc
  use variables, only : ngridx, ngridy, nx, ny, nlyr, &
       ntime, n_freq, kextcloud, KEXTATMO, &
       kextrain, kextice, kextsnow, kextgraupel, &
       salbtot, back, g_coeff, FREQ
  
  implicit none

  character(len=*), intent(in) :: OUT_FILE
  integer, intent(in) :: time_len, i_freq

  ! internal variables
  integer :: status, ncid, VarId
  integer :: nDims, unlimdimid
  character(len=len(OUT_FILE)+3) :: ncfile
  character(len=8) :: end_date
  character(len=10) :: end_time

  ncfile = OUT_FILE
  !status = NF90_OPEN_PAR(ncfile, IOR(IOR(NF90_WRITE, NF90_NETCDF4), &
  !     NF90_MPIIO), MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
  status = NF90_OPEN(ncfile, IOR(NF90_WRITE, NF90_NETCDF4), ncid)
  call check_nc(status, 'MP_STORENCDF() Opening the NetCDF', .true.)

  ! Getting NetCDF file dimensions and lengths:
  ! For frequency:
  if(i_freq.LT.1) stop 'ERROR: finding the frequency index in NP_STORENCDF'

  ! ---------------------------------------
  ! Writting the ATMOSPHERIC Extintion coeff
  status = nf90_inq_varid(ncid, "kext_atm", VarId)
  call check_nc(status, 'KEXT_ATMOS  ID cannot be read!', .true.)
  !status = nf90_var_par_access(ncid, VarId, NF90_INDEPENDENT)
  status = nf90_put_var(ncid, VarId, KEXTATMO, start=(/1, 1, 1, i_freq, time_len/) )
  call check_nc(status, 'KEXT_ATMOS cannot be written!', .true.)
  
  ! Writting the Cloud Extintion coeff
  status = nf90_inq_varid(ncid, "kext_qc", VarId)
  call check_nc(status, 'KEXT_CLOUD ID cannot be read!', .true.)
  !status = nf90_var_par_access(ncid, VarId, NF90_INDEPENDENT)
  status = nf90_put_var(ncid, VarId, kextcloud, start=(/1, 1, 1, i_freq, time_len/) )
  call check_nc(status, 'KEXT_CLOUD cannot be written!', .true.)

  ! Writting the Rain Extintion coeff
  status = nf90_inq_varid(ncid, "kext_qr", VarId)
  call check_nc(status, 'KEXT_RAIN ID cannot be read!', .true.)
  !status = nf90_var_par_access(ncid, VarId, NF90_INDEPENDENT)
  status = nf90_put_var(ncid, VarId, kextrain, start=(/1, 1, 1, i_freq, time_len/) )
  call check_nc(status, 'KEXT_RAIN cannot be written!', .true.)

  ! Writting the ICE Extintion coeff
  status = nf90_inq_varid(ncid, "kext_qi", VarId)
  call check_nc(status, 'KEXT_ICE ID cannot be read!', .true.)
  !status = nf90_var_par_access(ncid, VarId, NF90_INDEPENDENT)
  status = nf90_put_var(ncid, VarId, kextice, start=(/1, 1, 1, i_freq, time_len/) )
  call check_nc(status, 'KEXT_ICE cannot be written!', .true.)

  ! Writting the SNOW Extintion coeff
  status = nf90_inq_varid(ncid, "kext_qs", VarId)
  call check_nc(status, 'KEXT_SNOW ID cannot be read!', .true.)
  !status = nf90_var_par_access(ncid, VarId, NF90_INDEPENDENT)
  status = nf90_put_var(ncid, VarId, kextsnow, start=(/1, 1, 1, i_freq, time_len/) )
  call check_nc(status, 'KEXT_SNOW cannot be written!', .true.)

  ! Writting the GRAUPEL Extintion coeff
  status = nf90_inq_varid(ncid, "kext_qg", VarId)
  call check_nc(status, 'KEXT_GRAUPEL ID cannot be read!', .true.)
  !status = nf90_var_par_access(ncid, VarId, NF90_INDEPENDENT)
  status = nf90_put_var(ncid, VarId, kextgraupel, start=(/1, 1, 1, i_freq, time_len/) )
  call check_nc(status, 'KEXT_GRAUPEL cannot be written!', .true.)

  ! Writting the Total Albedo coeff
  status = nf90_inq_varid(ncid, "alb_tot", VarId)
  call check_nc(status, 'TOTAL Albedo ID cannot be read!', .true.)
  !status = nf90_var_par_access(ncid, VarId, NF90_INDEPENDENT)
  status = nf90_put_var(ncid, VarId, salbtot, start=(/1, 1, 1, i_freq, time_len/) )
  call check_nc(status, 'TOTAL  ALBEDO cannot be written!', .true.)

  ! Writting the Backscattering coeff
  status = nf90_inq_varid(ncid, "back_scatt", VarId)
  call check_nc(status, 'Backscatter ID cannot be read!', .true.)
  !status = nf90_var_par_access(ncid, VarId, NF90_INDEPENDENT)
  status = nf90_put_var(ncid, VarId, back, start=(/1, 1, 1, i_freq, time_len/) )
  call check_nc(status, 'Backscattering cannot be written!', .true.)

  ! Writting the asymetry factor
  status = nf90_inq_varid(ncid, "g_coeff", VarId)
  call check_nc(status, 'G-factor ID cannot be read!', .true.)
  !status = nf90_var_par_access(ncid, VarId, NF90_INDEPENDENT)
  status = nf90_put_var(ncid, VarId, g_coeff, start=(/1, 1, 1, i_freq, time_len/) )
  call check_nc(status, 'G-factor cannot be written!', .true.)

  if(time_len.EQ.NTIME.AND.i_freq.EQ.n_freq) then
     ! writting Initial date as global variable:
     call date_and_time(DATE=end_date, TIME=end_time)
     status = nf90_redef(ncid)
     status = nf90_put_att(ncid,NF90_GLOBAL,"End_Time", &
          end_date//' '//end_time )
     status = nf90_enddef(ncid)
  end if

  status = NF90_CLOSE(ncid)
  call check_nc(status, 'Closing NetCDF was not possible!', .true.)
  return
end subroutine MP_storencdf
! __________________________________________________________________/

! _____________________________________________________________________
! Subroutine to allocate/deallocate global variables
subroutine allocate_RT3_variables
  use variables
  ! ALLOCATING VARIABLES with netCDF dimensions:
  allocate( hgt_tmp(ngridx, ngridy, 0:nlyr, ntime) )
  allocate( temp_tmp(ngridx, ngridy, 0:nlyr, ntime) )
  allocate( press_tmp(ngridx, ngridy, 0:nlyr, ntime) )
  allocate( relhum_tmp(ngridx, ngridy, 0:nlyr, ntime) )
  allocate( mixr_tmp(ngridx, ngridy, nlyr, ntime) )
  allocate( cloud_water_tmp(ngridx, ngridy, nlyr, ntime) )
  allocate( rain_water_tmp(ngridx, ngridy, nlyr, ntime) )
  allocate( cloud_ice_tmp(ngridx, ngridy, nlyr, ntime) )
  allocate( snow_tmp(ngridx, ngridy, nlyr, ntime) )
  allocate( graupel_tmp(ngridx, ngridy, nlyr, ntime) )
  allocate( windvel_tmp(ngridx, ngridy, nlyr, ntime) )
  allocate( winddir_tmp(ngridx, ngridy, nlyr, ntime) )
  allocate( lat(ngridx, ngridy),  lon(ngridx, ngridy) )
  allocate( qidx(ngridx, ngridy, ntime) )
  allocate( UnixTime(ntime) )
  allocate( PBLH(ngridx, ngridy, ntime) )
  
  ! Initialazing the variables to a fixed value
  hgt_tmp = 0.0d0
  press_tmp = 0.0d0
  temp_tmp = 0.0d0
  relhum_tmp = 0.0d0
  mixr_tmp = 0.0d0
  cloud_water_tmp = 0.0d0
  rain_water_tmp = 0.0d0
  cloud_ice_tmp = 0.0d0
  snow_tmp = 0.0d0
  graupel_tmp = 0.0d0
  windvel_tmp = 0.0d0
  winddir_tmp = 0.0d0
  lat  = 0.0d0
  long = 0.0d0
  qidx = 0.0d0
  PBLH = -99.9
  
  return
end subroutine allocate_RT3_variables



! ---------------------------------------------------------------------
subroutine Broadcast_variables(error)
  use variables
  use mpi
  implicit none
      
  integer, intent(out) :: error
  integer :: rank, err=0

  call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)
 
  ! Broadcasting first the dimensions:
  call MPI_BCast(ngridx, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, err)
  call MPI_BCast(ngridy, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, err)
  call MPI_BCast(nlyr, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, err)
  call MPI_BCast(ntime, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, err) 
  if(rank.GT.0) then
     call allocate_RT3_variables
  end if

  ! After Allocating variables, broadcast them:
  call MPI_Bcast(hgt_tmp, size(hgt_tmp), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)
  call MPI_Bcast(temp_tmp, size(temp_tmp), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)
  call MPI_Bcast(press_tmp, size(press_tmp), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)
  call MPI_Bcast(relhum_tmp, size(relhum_tmp), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)
  call MPI_Bcast(mixr_tmp, size(mixr_tmp), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)
  call MPI_Bcast(cloud_water_tmp, size(cloud_water_tmp), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)
  call MPI_Bcast(rain_water_tmp, size(rain_water_tmp), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)
  call MPI_Bcast(cloud_ice_tmp, size(cloud_ice_tmp), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)
  call MPI_Bcast(snow_tmp, size(snow_tmp), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)
  call MPI_Bcast(graupel_tmp, size(graupel_tmp), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)
  call MPI_Bcast(windvel_tmp, size(windvel_tmp), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)
  call MPI_Bcast(winddir_tmp, size(winddir_tmp), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)
  call MPI_Bcast(lat, size(lat), MPI_REAL, 0, MPI_COMM_WORLD, err)
  call MPI_Bcast(lon, size(lon), MPI_REAL, 0, MPI_COMM_WORLD, err)
  call MPI_Bcast(qidx, size(qidx), MPI_INTEGER, 0, MPI_COMM_WORLD, err)
  call MPI_Bcast(UnixTime, size(UnixTime), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)
  call MPI_Bcast(PBLH, size(PBLH), MPI_REAL, 0, MPI_COMM_WORLD, err)

  error = err
  return
end subroutine Broadcast_variables




!!$subroutine MP_storencdf(OUT_FILE,time_len,i_freq,y_grid,x_grid,NLYR,LAYERS,TEMP,PRESS,RH,QV,QC,&
!!$     &WD, WS, KEXTQC, KEXTATM, KEXTTOT, ALBEDO, BACKSCATT, GCOEFF)
!!$
!!$  ! ********** Writting Station level Variables ****************************
!!$  ! Writting the 2m Temperature
!!$  status = nf90_inq_varid(ncid, "T2m", VarId)
!!$  if(status /= nf90_NoErr) stop 'T2m variable ID cannot be read!'
!!$  status = nf90_put_var(ncid, VarId, temp_lev(:,:,0:0), start=(/1 , 1, time_len/), count=(/ngridx, ngridy, 1/))
!!$  call check_nc(status, 'T2m variable cannot be written!', .true.)
!!$
!!$  ! Writting the 2m RH
!!$  status = nf90_inq_varid(ncid, "RH2m", VarId)
!!$  if(status /= nf90_NoErr) stop 'RH2m variable ID cannot be read!'
!!$  status = nf90_put_var(ncid, VarId, relhum_lev(:,:,0:0), start=(/1, 1, time_len/), count=(/ngridx, ngridy, 1/))
!!$  call check_nc(status, 'RH2m variable cannot be written!', .true.)
!!$  
!!$  ! Writting the 2m Air Pressure
!!$  status = nf90_inq_varid(ncid, "P2m", VarId)
!!$  if(status /= nf90_NoErr) stop 'Air pressure variable ID cannot be read!'
!!$  status = nf90_put_var(ncid, VarId, press_lev(:,:,0:0), start=(/1, 1, time_len/) )
!!$  call check_nc(status, 'P2m variable cannot be written!', .true.)
!!$  
!!$  ! ********** Writting Profile Variables ****************************
!!$  ! Writting the Temperature
!!$  status = nf90_inq_varid(ncid, "temp", VarId)
!!$  if(status /= nf90_NoErr) stop 'Temperature profile ID cannot be read!'
!!$  status = nf90_put_var(ncid, VarId, temp_lev(:,:,1:nlyr), start=(/1, 1, 1, time_len/) )
!!$  if(status /= nf90_NoErr) stop 'Temperature profile cannot be written!'
!!$  
!!$  ! Writting the Pressure
!!$  status = nf90_inq_varid(ncid, "press", VarId)
!!$  if(status /= nf90_NoErr) stop 'Pressure profile ID cannot be read!'
!!$  status = nf90_put_var(ncid, VarId, press_lev(:,:,1:nlyr), start=(/1, 1, 1, time_len/) )
!!$  if(status /= nf90_NoErr) stop 'Pressure profile cannot be written!'
!!$
!!$  ! Writting the Relative Humidity
!!$  status = nf90_inq_varid(ncid, "rh", VarId)
!!$  if(status /= nf90_NoErr) stop 'RH profile ID cannot be read!'
!!$  status = nf90_put_var(ncid, VarId, relhum_lev(:,:,1:nlyr), start=(/1, 1 ,1, time_len/) )
!!$  if(status /= nf90_NoErr) stop 'RH profile cannot be written!'
!!$
!!$  ! Writting the QV variable
!!$  status = nf90_inq_varid(ncid, "qv", VarId)
!!$  if(status /= nf90_NoErr) stop 'QV variable ID cannot be read!'
!!$  status = nf90_put_var(ncid, VarId, mixr_tmp(:,:,:,time_len), &
!!$       start=(/1, 1, 1, time_len/) )
!!$  if(status /= nf90_NoErr) stop 'QV cannot be written!'
!!$  ! Writting the QC variable
!!$  status = nf90_inq_varid(ncid, "qc", VarId)
!!$  if(status /= nf90_NoErr) stop 'QC variable ID cannot be read!'
!!$  status = nf90_put_var(ncid, VarId, cloud_water, start=(/1, 1, 1, time_len/) )
!!$  if(status /= nf90_NoErr) stop 'QC cannot be written!'
!!$
!!$  ! Writting the QR variable
!!$  status = nf90_inq_varid(ncid, "qr", VarId)
!!$  if(status /= nf90_NoErr) stop 'QR variable ID cannot be read!'
!!$  status = nf90_put_var(ncid, VarId, rain_water, start=(/1, 1, 1, time_len/) )
!!$  if(status /= nf90_NoErr) stop 'QR cannot be written!'
!!$
!!$  ! Writting the QI variable
!!$  status = nf90_inq_varid(ncid, "qi", VarId)
!!$  if(status /= nf90_NoErr) stop 'QI variable ID cannot be read!'
!!$  status = nf90_put_var(ncid, VarId, cloud_ice, start=(/1, 1, 1, time_len/) )
!!$  if(status /= nf90_NoErr) stop 'QI cannot be written!'
!!$
!!$  ! Writting the QS variable
!!$  status = nf90_inq_varid(ncid, "qs", VarId)
!!$  if(status /= nf90_NoErr) stop 'QS variable ID cannot be read!'
!!$  status = nf90_put_var(ncid, VarId, snow, start=(/1, 1, 1, time_len/) )
!!$  if(status /= nf90_NoErr) stop 'QS cannot be written!'
!!$
!!$  ! Writting the QG variable
!!$  status = nf90_inq_varid(ncid, "qg", VarId)
!!$  if(status /= nf90_NoErr) stop 'QG variable ID cannot be read!'
!!$  status = nf90_put_var(ncid, VarId, graupel, start=(/1, 1, 1, time_len/) )
!!$  if(status /= nf90_NoErr) stop 'QG cannot be written!'
!!$
!!$  ! Writting the Wind Direction variable
!!$  status = nf90_inq_varid(ncid, "wd", VarId)
!!$  if(status /= nf90_NoErr) stop 'WD variable ID cannot be read!'
!!$  status = nf90_put_var(ncid, VarId, winddir_tmp(:,:,:, time_len), start=(/1, 1, 1, time_len/) )
!!$  if(status /= nf90_NoErr) stop 'WD cannot be written!'
!!$  ! Writting the Wind Speed variable
!!$  status = nf90_inq_varid(ncid, "ws", VarId)
!!$  if(status /= nf90_NoErr) stop 'WS variable ID cannot be read!'
!!$  status = nf90_put_var(ncid, VarId, windvel_tmp(:,:,:, time_len), start=(/1, 1, 1, time_len/) )
!!$  if(status /= nf90_NoErr) stop 'WS cannot be written!'

