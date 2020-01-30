      module variables
      integer ngridx, ngridy, nlyr, ntime
      integer(kind=4), allocatable :: qidx(:,:,:)
      real(kind=8), allocatable, dimension(:,:,:,:) :: hgt_tmp,temp_tmp,
     $     press_tmp, relhum_tmp,
     $     mixr_tmp, cloud_water_tmp,
     $     rain_water_tmp, cloud_ice_tmp, snow_tmp, graupel_tmp,
     $     winddir_tmp, windvel_tmp
      real(kind=4) ,allocatable, dimension(:,:) :: lat, lon
      character(len=:), allocatable :: TimeStamp(:)

!     tests:
      real(kind=8), allocatable, dimension(:,:,:) :: kakes
!     Constants:
      real, parameter :: PI = dacos(-1.0d0)
      real, parameter :: PI2deg = 45.0/atan(1.0d0)
      end module variables
      
      program model_radtran_MW

C     Radiative transfer code to process COSMO-model, WRF-model derived profiles or
C     radiosonde measured atmospheric profiles.
C     The code read a full COSMO/WRF/RS grid and compute for each profile the 
c     radiative transfer for the given frequency
cC
C     This code is completely self contained and needs no databases
C     or lookup tables to run.  By convention, the quantities followed
C     by "_lev" are given at the layer heights while the qauntitites
c     w/o "_lev" are layer average quantities
C
C     NEW FEATURES: (added 2018 by Pablo Saavedra G.)
C     * forcing input as netCDF
C     * RT3 runs over the time dimension and grid dimension after one input read
C     * multiple frequency support
C     * parameters input as name-list (new) or standard input (old)
C     * one single output as netCDF file includes TBs and microphysics
C     * parallerism support via OpenMP (beta-version)
C     +-------+---------+---------+---------+---------+---------+---------+-+
C+----------------------------------------------------------------------
C                    **  RADTRAN I/O SPECIFICATIONS  **
      use variables
      implicit none
      include "omp_lib.h"    ! PSG: including OpenMP library
      include    'parameters.inc'
      integer i,j,k,isamp,jj,jsamp,nf,nx,ny,nz,
     $n_verify_input,nlev, !,ngridx,ngridy,nlyr
     $nx_in,nx_fin,ny_in,ny_fin,offset1,offset2,length1,length2
      character Nzstr*3,xstr*3,ystr*3,   ! PSG: Nzstr*2
     $ xstr1*3,ystr1*3,xstr2*3,ystr2*3,
     $ frq_str*5,theta_str*3 ,H_str*3,surf_type*10
      character input_file*99, micro_str*31,   ! PSG: original was input_file*99, frq_str*4
     $ SP_str*3,str1*1, DELTAM*1,
     $file_profile2*78,SD_snow*3,EM_snow*5,SD_grau*3,EM_grau*5,
     $ SD_rain*3,N0snowstr*3,N0graustr*3,N0rainstr*3
       REAL*8    DIRECT_FLUX, DIRECT_MU
      real*8      freq(20),lam,gammln     ! PSG: making FREQ vector for many frequencies
      character*2 year,month,day, hour    ! PSG: including hour
      character*8  date_str, sim_tag      ! PSG: was character*6  date_str, new sim_tag
      INTEGER AUIOF, BUIO                 ! PSG: open-file units depending on OMP_THREAD
       INTEGER     MAXLEG, NLEGEN, NUMRAD,NLEGENcw,NLEGENci,
     $NLEGENgr,NLEGENsn,NLEGENrr,aziorder, NUMAZIMUTHS,SRC_CODE
      REAL*8        RAD1, RAD2,refre,refim,SP
      REAL*8 N_0sr,Coeff_corr, AD,BD,ALPHA,GAMMA,emissivity, ! PSG: ,pi
     $ n0S,lambda_D,tmp,N_0snowD,N_0grauD,N_0snowDsnow,
     $ N_0grauDgrau,N_0rainD,Coeff_snow,a_mgraup, b_g, b_snow,
     $ a_msnow,
     $ Coeff_grau,LEGEN(200),LEGEN2(200),LEGEN3(200),LEGEN4(200),
     $LEGENcw(200),LEGENrr(200),LEGENci(200),
     $LEGENgr(200),LEGENsn(200),LEGENRain(200),
     $LEGEN2cw(200),LEGEN2rr(200),LEGEN2ci(200),
     $LEGEN2gr(200),LEGEN2sn(200), 
     $LEGEN3cw(200),LEGEN3rr(200),LEGEN3ci(200),
     $LEGEN3gr(200),LEGEN3sn(200), 
     $LEGEN4cw(200),LEGEN4rr(200),LEGEN4ci(200),
     $LEGEN4gr(200),LEGEN4sn(200),
     $  denliq,drop_mass,Deltar,denice,densnow,dengraup,
     $ fvol_ice
      COMPLEX*16    MINDEX,Im,m_air,m_MG,m_ice

      
      real*8       deltaxy(2),Tavg,ABSIND,ABSCOF  ! PSG: deltax, deltay -> deltaxy(2)
      real*8       tskin(mxgridx,mxgridy)
      integer    ics(mxgridx,mxgridy)
      real*8       sfcrain(mxgridx,mxgridy)
      real*8       hgt_lev(mxgridx,mxgridy,0:mxlyr)
      real*8       press_lev(mxgridx,mxgridy,0:mxlyr)
      real*8       temp_lev(mxgridx,mxgridy,0:mxlyr)
      real*8       relhum_lev(mxgridx,mxgridy,0:mxlyr)  
c      real*8       relhum(mxgridx,mxgridy,mxlyr)

      real*8       cloud_water(mxgridx,mxgridy,mxlyr)
      real*8       rain_water(mxgridx,mxgridy,mxlyr)
      real*8       max_rainwater(mxgridx,mxgridy)
      real*8       cloud_ice(mxgridx,mxgridy,mxlyr)
      real*8       snow(mxgridx,mxgridy,mxlyr)
      real*8       graupel(mxgridx,mxgridy,mxlyr)
  
      real*8       cloud_water_lev(mxgridx,mxgridy,0:mxlyr)
      real*8       rain_water_lev(mxgridx,mxgridy,0:mxlyr)
      real*8       cloud_ice_lev(mxgridx,mxgridy,0:mxlyr)
      real*8       snow_lev(mxgridx,mxgridy,0:mxlyr)
      real*8       graupel_lev(mxgridx,mxgridy,0:mxlyr)
    
      real*8        LYR_TEMP(0:MXLYR), LYR_PRES(0:MXLYR),
     $REL_HUM(MXLYR),KEXTATMO(MXLYR),AVGPRESSURE(MXLYR),
     $ VAPORPRESSURE(MXLYR),P11(2),ang(2)


      real*8       atm_ext, kextcw, salbcw, asymcw, kextrr, salbrr,
     $           asymrr, kextci, salbci, asymci, kextsn, salbsn, 
     $           asymsn, kextgr, salbgr, asymgr, salbhl,
     $           asymhl,   backcw, backrr, 
     $          backci, backsn, backgr,tau_min,tau_max
      real*8       avg_pressure(mxgridx,mxgridy,mxlyr),
     $           vapor_pressure(mxgridx,mxgridy,mxlyr),
     $         rho_vap(mxgridx,mxgridy,mxlyr)
      real*8       kexttot(mxgridx,mxgridy,mxlyr),   
     $ g_coeff(mxgridx,mxgridy,mxlyr),
     $ kextcloud(mxgridx,mxgridy,mxlyr),
     $ kextrain(mxgridx,mxgridy,mxlyr),  
     $ kextice(mxgridx,mxgridy,mxlyr),   
     $ kextgraupel(mxgridx,mxgridy,mxlyr),   
     $ kextsnow(mxgridx,mxgridy,mxlyr),
     $  salbtot(mxgridx,mxgridy,mxlyr), 
     $  tau(mxgridx,mxgridy), tau_hydro(mxgridx,mxgridy),
     $  abs(mxgridx,mxgridy,mxlyr),
     $           asymtot(mxgridx,mxgridy,mxlyr),
     $              back(mxgridx,mxgridy,mxlyr),
     $mu,D_0,D_max,D_min,A1,A2
      integer j_temp,ind_temp,N_lay_cut
 
      INTEGER  MAXV,MAXA,MAXLAY !,NUMMU,ntheta_i,ntheta_s
      PARAMETER (MAXV=90, MAXA=32) ! PSG: (MAXV=64, MAXA=32)
      PARAMETER (MAXLAY=200)
       INTEGER  NSTOKES
      INTEGER  NOUTLEVELS, OUTLEVELS(MAXLAY)
      CHARACTER QUAD_TYPE*1, UNITS*1, OUTPOL*2, GROUND_TYPE*1,
     $rLWC_str*4
    
       CHARACTER*164 OUT_FILE,file_PH(mxlyr), NCDFOUT,
     $ file_PH2(mxlyr),tmp_file1, origin_str   ! PSG: change *64 to *164 and added NCDFOUT & origin_str
       CHARACTER*68 file_profile
      REAL*8   MU_VALUES(MAXV), GROUND_TEMP, GROUND_ALBEDO
      REAL*8   SKY_TEMP, WAVELENGTH
       COMPLEX*16  GROUND_INDEX

       real*8  TTheta_i(ntheta_i),TTheta_s(ntheta_i),Salinity

      
         character ssstr*1,ttstr*1,Anglestr*4,FILEOUT3D*65
      integer i_bot,i_top,nnz,N_layer_new,ss,tt,
     $ length
       real*8 H_levs(mxlyr+1),Angle_view,Angle_zenith,
     $I1,I2,Angle_viewdeg,Upar,Vpar,
     $Ds_cloud_slant,Ds_obs_point,DH_bot_intersect,H_bot,
     $H_top,DH_top_intersect,K_extbot,K_exttop,T_bot,T_top
 
C     !PSG: Following temporal varaibles is to include NetCDF time depending
       integer n_freq, ifreq, timeidx, i_time    ! PSG: new block... , ntime
       character frq_idx*5
       !real*8 hgt_tmp(mxgridx,mxgridy,0:mxlyr,mxtime)
       !real*8 press_tmp(mxgridx,mxgridy,0:mxlyr,mxtime)
       !real*8 temp_tmp(mxgridx,mxgridy,0:mxlyr,mxtime) ! PSG: temp for module
       !real*8 relhum_tmp(mxgridx,mxgridy,0:mxlyr,mxtime)
       !real*8 mixr_tmp(mxgridx,mxgridy,0:mxlyr,mxtime)
       !real*8 cloud_water_tmp(mxgridx,mxgridy,mxlyr,mxtime)
       !real*8 rain_water_tmp(mxgridx,mxgridy,mxlyr,mxtime)
       !real*8 cloud_ice_tmp(mxgridx,mxgridy,mxlyr,mxtime)
       !real*8 snow_tmp(mxgridx,mxgridy,mxlyr,mxtime)
       !real*8 graupel_tmp(mxgridx,mxgridy,mxlyr,mxtime)
       !real*8 winddir_tmp(mxgridx,mxgridy,mxlyr,mxtime)
       !real*8 windvel_tmp(mxgridx,mxgridy,mxlyr,mxtime)
       !integer*4 qidx(mxgridx,mxgridy,mxtime)
       real*4 yy(mxgridx,mxgridy,mxtime), mm(mxgridx,mxgridy,mxtime)
       real*4 dd(mxgridx,mxgridy,mxtime), hh(mxgridx,mxgridy,mxtime)
       !real*4 lat(mxgridx,mxgridy),lon(mxgridx,mxgridy)

       namelist/inputparam/input_file, nx_in,nx_fin,ny_in,ny_fin,n_freq,
     $      tau_min,tau_max, FREQ,
     $      GAS_EXTINCTION, SD_snow, N_0snowDsnow,EM_snow,SP,SD_grau,
     $      N_0grauDgrau, EM_grau,SD_rain, N_0rainD
       integer istatus

C     !PSG: -- end of definition for NetCDF temporal variables

      LOGICAL     PHASEFLAG
   
        LOGICAL GAS_EXTINCTION
         REAL*8 E1,E2
       CHARACTER*100 FILEOUT1
c     The following parameters are set for the TMI sensor

       Im=(0.0d0,1.0d0)
       !pi = dacos(-1.0d0)
c    some inputs variable  
       Nstokes=2
C       N_lay_cut=135  ! PSG: commented, not needed
       QUAD_TYPE='D' ! PSG: original was 'D'
       GROUND_TYPE='F'  !'L'
       Aziorder=0
       NUMAZIMUTHS=1
        DELTAM='N'
       SRC_CODE=2
       DIRECT_FLUX=0.d0
       DIRECT_MU=0.0d0       
       maxleg=200
        emissivity=0.90d0
        GROUND_ALBEDO=1.0-emissivity
c   
       SKY_TEMP=2.7
C       WAVELENGTH=1000000.0    ! PSG: must be freq dependent :/
       UNITS='T'
       OUTPOL='VH'
       NOUTLEVELS=2
c     write(11,*),'GAS_EXTINCT',OUTLEVELS(1),OUTLEVELS(2)
c      stop
c 
c     Get input/output file names from command line.  
c
!     PSG: input block for input parameters:
       OPEN(UNIT=100,FILE='input',STATUS='old',IOSTAT=istatus)
       if(istatus.eq.0) then
          READ(UNIT=100,nml=inputparam)
          close(UNIT=100)
       else
          WRITE (*,'(1X,A)') 'CRM input file'
          READ (*,*) input_file
c      input_file='/work/batta/colo/CRM/TCOF22_GCE01.00480.profile'
c       input_file='/work/batta/colo/CRM/MIDACF_GCE01.00604.profile'
c      input_file='/work/batta/colo/CRM/LBAEST_GCE01.00360.profile'
          WRITE (*,'(1X,A)') 'nx_in,ny_in,ny_in,ny_fin'
          read(*,*) nx_in,nx_fin,ny_in,ny_fin,n_freq,tau_min,tau_max
c     WRITE (18,*) nx_in,ny_in,ny_in,ny_fin,tau_min,tau_max      
          WRITE (*,'(1X,A)') 'WHICH frequency (GHz)?'
          READ (*,*) FREQ(1:n_freq)
   
          WRITE (*,'(1X,A)')  'DO YOU WANT TO COMPUTE THE EXTINCTION
     $FROM ATMOSPHERIC GASES AS WELL?'
          READ (*,*)  GAS_EXTINCTION

          READ (*,*) SD_snow, N_0snowDsnow,EM_snow,SP,SD_grau,
     $         N_0grauDgrau,
     $         EM_grau,SD_rain, N_0rainD !N_0 are for the distr of snow/grau diameter
       endif   ! PSG: end for input parameters

        PHASEFLAG=.true.
C        LAM=299.7925/freq !mm 
                
        print*,'netCDF input files is '//input_file
C     !PSG: Calling the NetCDF routine to read data
        call read_wrf(len_trim(input_file), input_file,
     $       deltaxy, origin_str)

        !print*, 'mysu= ', ngridx, ngridy, nlyr, ntime

        !do i=1, ntime
        !   write(*,'(10F9.3)') (cloud_ice_tmp(i,j,5,1), j=1,mxgridy)
        !   print*, TimeStamp(i)
        !   write(*,'(10F9.3)') (temp_tmp(i,j,5,1), j=1,mxgridy)
        !enddo
        !print*, 'press at ground is: ', origin_str

        OUTLEVELS(1)=1          ! PSG: moved from befor call RT3 to here
        OUTLEVELS(2)=nlyr+1     ! PSG: N_lay_cut+1

        
        write(SP_str(1:3),'(f3.1)')SP
c     write(18,*)'str',H_str

        if (N_0snowDsnow.le.9.95d0) then
           write(N0snowstr,'(f3.1)')N_0snowDsnow
        else
           write(N0snowstr,'(f3.0)')N_0snowDsnow
        endif 
        
        if (N_0rainD.le.9.95) then
           write(N0rainstr,'(f3.1)')N_0rainD
        else
           write(N0rainstr,'(f3.0)')N_0rainD
        endif      
        if (N_0grauDgrau.le.9.95) then
           write(N0graustr,'(f3.1)')N_0grauDgrau
        else
           write(N0graustr,'(f3.0)')N_0grauDgrau
        endif               

        micro_str=SD_snow//N0snowstr//EM_snow//SP_str//
     $       SD_grau//N0graustr//EM_grau//SD_rain//N0rainstr

C     NetCDF output file definition (creating one Ncdf file per frequency):
        j_temp = scan(input_file,"/",back=.true.)+1
        
        NCDFOUT='../output/TB/RT3TB_'//trim(input_file(j_temp:)) !micro_str//'_'

        write(*,*) 'Creating NetCDF '//trim(NCDFOUT)
        write(*,*) 'ntime=',ntime,'; nlayer=',nlyr,'; nfreq=',n_freq

        call createncdf(len_trim(NCDFOUT),trim(NCDFOUT),
     $       NUMMU,n_freq,NSTOKES,nlyr,
     $       ngridx,ngridy,hgt_tmp(1,1,1:nlyr,1),
     $       FREQ(1:n_freq),input_file,micro_str,
     $       real(hgt_tmp(:ngridx,:ngridy,0,1),4),
     $       lat,lon,trim(origin_str))

C     !PSG: -- end of NetCDF reading routine
        do 777 ifreq=1,n_freq   ! PSG: include Frequency loop
           if (freq(ifreq).gt.100.d0) then
              write(frq_str,'(f6.2)') FREQ(ifreq)
              write(frq_idx,'(I6.6)') ifreq
           else
              write(frq_str,'(f5.2)') FREQ(ifreq)
              write(frq_idx,'(I5.5)') ifreq
           endif

c     write(*,29) frq_str  
           
           call omp_set_num_threads(4)
!$OMP PARALLEL NUM_THREADS(1) PRIVAD(AUIOF,BUIOF)
!$OMP DO
        do 656 nx=1, 2 !ngridx ! nx_in,nx_fin
          write(xstr,'(i3.3)') nx
                   
         do 656 ny = 1, 2 !ngridy  ! ny_in,ny_fin
            write(ystr,'(i3.3)') ny

C     !PSG: Passing temporal variables to old variables (no time)
            i_time = 0

            DO 646, timeidx = 1, 1 !ntime
               write(*,*) 'running on thread: ', OMP_GET_THREAD_NUM(),
     $              OMP_GET_MAX_THREADS()

!     Tests:
               allocate(kakes(ngridx, ngridy, 0:nlyr))
               kakes = temp_tmp(:,:,:,timeidx)
               print*, shape(kakes)
               print*, shape(temp_tmp)
               write(*,'(11F9.2)') (kakes(31,1,k), k=0,10)
               write(*,'(11F9.2)') (temp_tmp(31,1,k,timeidx), k=0,10)
               stop
               hgt_lev = hgt_tmp(:,:,:, timeidx)
           press_lev = press_tmp(:,:,:, timeidx)  ! PSG: i_time)
           temp_lev = temp_tmp(:,:,:, timeidx)
           relhum_lev = relhum_tmp(:,:,:, timeidx)
           cloud_water = cloud_water_tmp(:,:,:, timeidx)
           rain_water = rain_water_tmp(:,:,:, timeidx)
           cloud_ice = cloud_ice_tmp(:,:,:, timeidx)
           snow = snow_tmp(:,:,:, timeidx)
           graupel = graupel_tmp(:,:,:, timeidx)
           max_rainwater = maxval(rain_water,DIM=3)

c     
c     Read the standard ASCII COSMO cloud  model output
c
c$$$
c$$$      OPEN(UNIT = 14, FILE = input_file, STATUS = 'OLD',
c$$$     $     form = 'formatted')
c$$$     
c$$$      read(14,*) year,month,day,ngridx,ngridy,nlyr,deltaxy(1),deltaxy(2)
c$$$C     write(18,*) year,month,day,ngridx,ngridy,nlyr,deltaxy(1),deltaxy(2)
c$$$      do i=1,2 !ngridx
c$$$         do j=1,2 !ngridy
c$$$            max_rainwater(i,j)=0.d0
c$$$            read(14,*) isamp, jsamp
c$$$            write(28,*) isamp, jsamp
c$$$            read(14,*) hgt_lev(i,j,0), press_lev(i,j,0),
c$$$     +           temp_lev(i,j,0), relhum_lev(i,j,0)
c$$$C            press_lev(i,j,0)=press_lev(i,j,0)/100.d0   ! PSG: input already in [hPa]
c$$$            do k=1,nlyr
c$$$               read(14,*) hgt_lev(i,j,k), press_lev(i,j,k),
c$$$     +              temp_lev(i,j,k), relhum_lev(i,j,k),
c$$$     +              cloud_water(i,j,k),rain_water(i,j,k),
c$$$     +              cloud_ice(i,j,k),snow(i,j,k), graupel(i,j,k)
c$$$               max_rainwater(i,j)=max(max_rainwater(i,j),
c$$$     +              rain_water(i,j,k))
c$$$C               press_lev(i,j,k)=press_lev(i,j,k)/100.d0   ! PSG: input already in [hPa]
c$$$            end do
C            write(*,*) "showing the readings..."
C            write(*,*) hgt_lev(i,j,0), press_lev(i,j,0),
C     +           temp_lev(i,j,0), relhum_lev(i,j,0), nlyr
C            do k=1,nlyr
C               write(*,*) k,hgt_lev(i,j,k), press_lev(i,j,k),
C     +              temp_lev(i,j,k), relhum_lev(i,j,k)
C            end do

c     do k=1,nlyr
c        cloud_water(i,j,k)=0.5*(cloud_water_lev(i,j,k)+ 
c     $cloud_water_lev(i,j,k-1))
c        rain_water(i,j,k)=0.5*(rain_water_lev(i,j,k)+ 
c     $rain_water_lev(i,j,k-1))
c       cloud_ice(i,j,k)=0.5*(cloud_ice_lev(i,j,k)+ 
c     $cloud_ice_lev(i,j,k-1))
c       snow(i,j,k)=0.5*(snow_lev(i,j,k)+ 
c     $snow_lev(i,j,k-1))
c       graupel(i,j,k)=0.5*(graupel_lev(i,j,k)+ 
c     $graupel_lev(i,j,k-1))
c          end do
c$$$        end do
c$$$      end do

c   503 format(f5.1,3f7.2,6f7.4,f7.3,5f7.2)
      
c
c     This GCE model format does not have all the fields expected by
c     the radiative transfer code (i.e. total pressure, and water vapor
c     pressure for this model).  Assign/compute the missing fields first. 
c
      call get_atmosG0(temp_lev,press_lev,relhum_lev,ngridx,
     +   ngridy,nlyr,avg_pressure, vapor_pressure,rho_vap)
      
c$$$      do 656 ifreq=1,n_freq                  ! PSG: include Frequency loop
c$$$        if (freq(ifreq).gt.100.d0) then
c$$$        write(frq_str,'(f5.1)') FREQ(ifreq)
c$$$        else
c$$$        write(frq_str,'(f4.1)') FREQ(ifreq)
c$$$        endif
c$$$c        write(*,29) frq_str  
c$$$          write(SP_str(1:3),'(f3.1)')SP
c$$$c         write(18,*)'str',H_str

c
c       Calculate optical properties of each volume element of
c       the cloud.  Optionally, write this file out for future use
c       by this code or another code.
c

          LAM=299.7925/freq(ifreq) !mm  :PSG
          WAVELENGTH = LAM*1E3    !micrometre  :PSG
          alpha=0.0d0
          gamma=1.0d0 ! always exponential SD 
          
          write(xstr1,'(i3.3)') nx_in 
          write(xstr2,'(i3.3)') nx_fin
          write(ystr1,'(i3.3)') ny_in
          write(ystr2,'(i3.3)') ny_fin
c            if (day.le.9) then
c          write(daystr,'(i1)') day
c          daystr='0'//daystr
c          else
c          write(daystr,'(i2)')day
c          endif 
c          write(18,*)'day',day,daystr
c          stop         


c$$$          if (N_0snowDsnow.le.9.95d0) then
c$$$          write(N0snowstr,'(f3.1)')N_0snowDsnow
c$$$          else
c$$$          write(N0snowstr,'(f3.0)')N_0snowDsnow
c$$$          endif 
c$$$ 
c$$$          if (N_0rainD.le.9.95) then
c$$$          write(N0rainstr,'(f3.1)')N_0rainD
c$$$          else
c$$$          write(N0rainstr,'(f3.0)')N_0rainD
c$$$          endif      
c$$$           if (N_0grauDgrau.le.9.95) then
c$$$          write(N0graustr,'(f3.1)')N_0grauDgrau
c$$$          else
c$$$          write(N0graustr,'(f3.0)')N_0grauDgrau
c$$$          endif               
c$$$
c$$$          micro_str=SD_snow//N0snowstr//EM_snow//SP_str//
c$$$     $ SD_grau//N0graustr//EM_grau//SD_rain//N0rainstr
      
          
        
            ! Checking whether the profile passed the qualitity control:
          if(qidx(nx, ny, timeidx).NE.15) go to 646

            i_time = i_time + 1
            !write(*,*) 'running on time: ',i_time
            !!write(year,'(I2.2)') int(yy(nx,ny, timeidx)-2000)
            !!write(month,'(I2.2)') int(mm(nx,ny, timeidx))
            !!write(day,'(I2.2)') int(dd(nx,ny, timeidx))
            !!write(hour,'(I2.2)') int(hh(nx,ny, timeidx))
            year = TimeStamp(timeidx)(3:4)
            month = TimeStamp(timeidx)(6:7)
            day = TimeStamp(timeidx)(9:10)
            hour = TimeStamp(timeidx)(12:13)
            
            date_str=year//month//day//hour 
c   computing the refractive index of the sea surface
        Salinity=10.0d0
        call DIECON(Salinity,temp_lev(nx,ny,0)-273.16,FREQ(ifreq),E1,E2)
        GROUND_INDEX=dconjg(sqrt(E1+Im*E2))
C     write(18,*)'maxrain',GROUND_INDEX,E1,E2

c             if (max_rainwater(nx,ny).le.1e-1) goto 646
             tau(nx,ny)=0.0d0 
             tau_hydro(nx,ny)=0.0d0
              do i=1,NLYR
             LYR_TEMP(0)=temp_lev(nx,ny,0)
              LYR_PRES(0)=0.1*press_lev(nx,ny,0) 
              LYR_TEMP(i)=temp_lev(nx,ny,i)
              LYR_PRES(i)=0.1*press_lev(nx,ny,i) 
              REL_HUM(i)= 0.5*(relhum_lev(nx,ny,i)+
     $ relhum_lev(nx,ny,i-1))
              enddo
              GROUND_TEMP= LYR_TEMP(0)
            IF (GAS_EXTINCTION) THEN
      CALL GET_ATMOSG(LYR_TEMP,LYR_PRES, REL_HUM, 
     $ MXLYR, NLYR, AVGPRESSURE, VAPORPRESSURE,FREQ(ifreq),KEXTATMO)
       ELSE 
         DO 30 I = 1,NLYR
          KEXTATMO(I)=0.0D0
 30          CONTINUE
       ENDIF   

C     PSG: add 'tmp/'
       file_profile='../output/tmp/Profilex'
       file_profile=trim(file_profile)//xstr//'y'//ystr//'f'//frq_str
C     file_profile='Profilex'//xstr//'y'//ystr//'f'//
C     $ frq_str 
    
        
        

c 1017          format(i2,1x,f3.1,1x,f5.1,1x,f6.1)
c          
cc     INITIALIZATION OF LEGENDRE COEFFICIENTS  CCCCCCCCCCC          
       do 1009 nz = 1,nlyr      ! PSG:  N_lay_cut  
          write(Nzstr,'(i3.3)') Nz ! PSG: '(i2.2)'

             NLegen=0           
             NLEGENcw=0
             NLEGENci=0
             NLEGENrr=0
             NLEGENsn=0
             NLEGENgr=0
             FILE_PH(nz) =''      !strings with blank spaces
            
c              write(18,*) 'nz',nz
c            write(2,1002) nz, hgt_lev(nx,ny,nz),temp_lev(nx,ny,nz),
c     $ 0.1*press_lev(nx,ny,nz), relhum(nx,ny,nz)
c 1002    format(i2,1x,f5.1,1x,f5.1,1x,f6.1,1x,f5.1)   
                do jj=1,200
              LEGEN(jj)=0.0d0
              LEGENcw(jj)=0.0d0
              LEGENci(jj)=0.0d0
              LEGENrr(jj)=0.0d0
              LEGENgr(jj)=0.0d0
              LEGENsn(jj)=0.0d0
        

              LEGEN2(jj)=0.0d0
              LEGEN2cw(jj)=0.0d0
              LEGEN2ci(jj)=0.0d0
              LEGEN2rr(jj)=0.0d0
              LEGEN2gr(jj)=0.0d0
              LEGEN2sn(jj)=0.0d0
          

              LEGEN3(jj)=0.0d0
              LEGEN3cw(jj)=0.0d0
              LEGEN3ci(jj)=0.0d0
              LEGEN3rr(jj)=0.0d0
              LEGEN3gr(jj)=0.0d0
              LEGEN3sn(jj)=0.0d0
         
      
              LEGEN4(jj)=0.0d0
              LEGEN4cw(jj)=0.0d0
              LEGEN4ci(jj)=0.0d0
              LEGEN4rr(jj)=0.0d0
              LEGEN4gr(jj)=0.0d0
              LEGEN4sn(jj)=0.0d0
          
              enddo

c     computing mean T of the layer
            Tavg = 0.5*(temp_lev(nx,ny,nz-1)+temp_lev(nx,ny,nz))

CCCCCCCCCCCCCC    SINGLE SCATTERING PROPERTIES OF CLOUD WATER  CCCCCCCCCCCCCCCCCCCCCCCC
c        write(*,*)'entro clw',freq(nf),Tavg,cloud_water(nx,ny,nz)
         if(cloud_water(nx,ny,nz).ge.1e-5) then  

         call REFWAT(1,LAM,Tavg,Refre,refim,ABSIND,ABSCOF) 
c          write(18,*)'ref',Refre,refim 
         mindex=refre-Im*refim
          Deltar=0.00001d0  !mm
          rad1=0.01d0     !10 micron radius monodisperse
          rad2=rad1+Deltar
          denliq=1.0d0
          drop_mass =4./3.*pi*rad1*rad1*rad1*denliq
          AD=cloud_water(nx,ny,nz)/drop_mass/Deltar
          BD=0.0d0
          alpha=0.0d0
          gamma=1.0d0 ! exponential SD 

          numrad=2
         
            call MIE (LAM, MINDEX, RAD1, RAD2, NUMRAD, MAXLEG,
     .                AD, BD, ALPHA, GAMMA, PHASEFLAG,
     .                 kextcw, salbcw,backcw, NLEGENcw, LEGENcw,
     $ LEGEN2cw,LEGEN3cw,LEGEN4cw, 'G')
           Nlegen=max(Nlegen,NLEGENcw)
c            write(18,*)'esco cw2', kextcw, salbcw, legen(1),legen(2),
c     $legen(3),legen(4),legen(20),mindex
          else
           kextcw=0.0d0
           salbcw=0.0d0
           backcw=0.0d0
          endif
CCCCCCCCCCCCCC    SINGLE SCATTERING PROPERTIES OF RAIN  CCCCCCCCCCCCCCCCCCCCCCCC

          if(rain_water(nx,ny,nz).gt.1e-5) then
c            call mie_rain(FREQ, Tavg, rain_water(nx,ny,nz),
c     $                   kextrr, salbrr, asymrr)
c            write(18,*)'esco rainm', kextrr, salbrr, asymrr
c          call refr_indx_liq(FREQ, Tavg, 1.0d0, refre, refim ) 
c          write(18,*)'ref',Refre,refim 
          
        call REFWAT(1,LAM,Tavg,Refre,refim,ABSIND,ABSCOF)
c          write(18,*)'ref',Refre,refim
         mindex=refre-Im*refim
          denliq=1.0d0
          rad1=0.005d0  !minimum radius in mm
          rad2=4.50d0   !maximum radius
          AD=2.0d0*N_0rainD    !this are for integration over radii (not diameters)!
          BD=2.d0*(pi*denliq*N_0rainD/rain_water(nx,ny,nz))**0.25
        
          numrad=100
        
            call MIE (LAM, MINDEX, RAD1, RAD2, NUMRAD, MAXLEG,
     .                AD, BD, ALPHA, GAMMA, PHASEFLAG,
     .                 kextrr, salbrr, backrr,NLEGENrr, LEGENrr,
     $ LEGEN2rr,LEGEN3rr,LEGEN4rr, 'G')
        Nlegen=max(Nlegen,NLEGENrr)
        
           else
           kextrr=0.0d0
           salbrr=0.0d0 
           backrr=0.0d0
            endif
C$$$$$$$$$$$$$$$$$ SINGLE SCATTERING PROPERTIES OF ICE CRYSTALS   $$$$$$$$$$$$$$$$$$$$$$$$$$$

            if(cloud_ice(nx,ny,nz).ge.1e-5) then
c            call mie_ci( FREQ, Tavg, cloud_ice(nx,ny,nz),
c     $                   kextci, salbci, asymci)
c               write(18,*)'esco ci', kextci, salbci, asymci
c             call refr_indx_ice(FREQ, Tavg, refre, refim)
c              write(18,*)'ref ice',Refre,refim 
             call REFICE(1,LAM,Tavg,Refre,refim,ABSIND,ABSCOF) 
c              write(18,*)'ref',Refre,refim 
              mindex=refre-Im*refim
              Deltar=0.0001d0   !mimicking a monodisperse distribution
              rad1=0.01d0    !10 micron radius
              rad2=rad1+Deltar
              denice=0.917d0
              drop_mass =4./3.*pi*rad1*rad1*rad1*denice
              AD=cloud_ice(nx,ny,nz)/drop_mass/Deltar
              BD=0.0d0
              alpha=0.0d0
              gamma=1.0d0 ! exponential SD 
              numrad=2
              call MIE (LAM, MINDEX, RAD1, RAD2, NUMRAD, MAXLEG,
     .                AD, BD, ALPHA, GAMMA, PHASEFLAG,
     .                 kextci, salbci, backci,NLEGENci, LEGENci, 
     $ LEGEN2ci, LEGEN3ci, LEGEN4ci, 'G')
              Nlegen=max(Nlegen,NLEGENci)
c          write(18,*)'esco ci2',kextci,salbci,legen(1),legen(2),
c     $legen(3),legen(4),legen(20),Nlegen,mindex
             else
           kextci=0.0d0
           salbci=0.0d0 
           backci=0.0d0
             endif

C$$$$$$$$$$$$$$$$$ SINGLE SCATTERING PROPERTIES OF SNOW  CCCCCCCCCCCCCCCCCCCCC
            if(snow(nx,ny,nz).ge.1e-5) then
            b_snow=2.0d0    !MKS system 
            a_msnow=0.038d0   

        
         call REFICE(1,LAM,Tavg,Refre,refim,ABSIND,ABSCOF) 
           m_ice=refre+Im*refim
           m_air=1.0d0+0.0d0*Im
          rad1=0.01
          !DSD are expressed in terms of snow radii 
         AD=2.0d0*N_0snowDsnow*dexp(0.107*(273.15-Tavg)) !Field param. ! multiplied by 10^6 is 1/m^4
         BD=2e-3*(dexp(Gammln(b_snow+1))*a_msnow*N_0snowDsnow
     $ *dexp(0.107*(273.15-Tavg))*1e6/
     $ (1e-3*snow(nx,ny,nz)))**(1.0d0/(1.0d0+b_snow))  !mm^-1   !formula 3.12 Mario Mech´s but for radii and units converted
         rad2=15.0/BD     
          numrad=500

         IF(EM_snow.eq.'softs') THEN
           call MIE_densitysizedep_softsphere(LAM,M_ice,m_air,
     $  a_msnow, b_snow, RAD1, RAD2, NUMRAD, MAXLEG,
     .                AD, BD, ALPHA, GAMMA, PHASEFLAG,
     .                  kextsn, salbsn,backsn, NLEGENsn, LEGENsn,
     $ LEGEN2sn,LEGEN3sn,LEGEN4sn,'G')
c       write(18,*)'new values snow', kextsn, salbsn, backsn

        ELSEIF(EM_snow.eq.'icesf') THEN 
          call MIE_densitysizedep_spheremasseq(LAM,M_ice,m_air,
     $  a_msnow, b_snow, RAD1, RAD2, NUMRAD, MAXLEG,
     .                AD, BD, ALPHA, GAMMA, PHASEFLAG,
     .                  kextsn, salbsn,backsn, NLEGENsn, LEGENsn,
     $ LEGEN2sn,LEGEN3sn,LEGEN4sn,'G') 
         ELSEIF(EM_snow.eq.'snowA'.or.EM_snow.eq.'snowB'.or.
     $EM_snow.eq.'roset'.or.EM_snow(1:3).eq.'Kim' ) THEN 
      call MIE_densitysizedep_parameterization(LAM,M_ice,m_air,
     $a_msnow,b_snow,EM_snow,RAD1,RAD2,NUMRAD,MAXLEG,
     .AD, BD, ALPHA, GAMMA, PHASEFLAG,
     .                  kextsn, salbsn,backsn, NLEGENsn, LEGENsn,
     $ LEGEN2sn,LEGEN3sn,LEGEN4sn,'G')

           ELSE
           write(18,*)'no em mod',EM_snow
           stop 
           ENDIF
            Nlegen=max(Nlegen,NLEGENsn)
       call LEGEndre2PHASEFUNCTION(LEGENsn,NLEGENsn,
     $ 2,200,P11,ANG) 
            backsn=kextsn*salbsn*P11(2)     
c          write(18,*)'esco sn2',snow(nx,ny,nz),densnow,
c     $ kextsn,salbsn,legensn(1),legensn(2),
c     $legensn(3),legensn(4),legensn(20),Nlegen,mindex,nx
              else
 444           continue
            kextsn=0.0d0
            salbsn=0.0d0
            backsn=0.0d0 
             endif
          
C$$$$$$$$$$$$$$$$$ SINGLE SCATTERING PROPERTIES OF GRAUPEL CCCCCCCCCCCCCCCCCCCC
             if(graupel(nx,ny,nz).ge.1e-5) then
c             dengraup=0.4
            b_g=3.1d0
            a_mgraup=169.6d0



         
         call REFICE(1,LAM,Tavg,Refre,refim,ABSIND,ABSCOF) 
           m_ice=refre+Im*refim
           m_air=1.0d0+0.0d0*Im
          rad1=0.01
          numrad=500

         AD=2.0d0*N_0grauDgrau
         BD=2e-3*(dexp(Gammln(b_g+1))*a_mgraup*N_0grauDgrau*1e6/
     $ (1e-3*graupel(nx,ny,nz)))**(1.0d0/(1.0d0+b_g))  !mm^-1
         rad2=15.0/BD

           IF(EM_grau.eq.'softs') THEN
         call MIE_densitysizedep_softsphere (LAM,M_ice,m_air,
     $  a_mgraup, b_g, RAD1, RAD2, NUMRAD, MAXLEG,
     .                AD, BD, ALPHA, GAMMA, PHASEFLAG,
     .                 kextgr, salbgr, backgr,NLEGENgr, LEGENgr, 
     $ LEGEN2gr, LEGEN3gr, LEGEN4gr, 'G')

            ELSEIF(EM_grau.eq.'icesf') THEN 
           call MIE_densitysizedep_spheremasseq(LAM, M_ice,m_air,
     $  a_mgraup, b_g, RAD1, RAD2, NUMRAD, MAXLEG,
     .                AD, BD, ALPHA, GAMMA, PHASEFLAG,
     .                 kextgr, salbgr, backgr,NLEGENgr, LEGENgr, 
     $ LEGEN2gr, LEGEN3gr, LEGEN4gr, 'G')

c       write(18,*)'new values', kextgr, salbgr, backgr


          ELSE
           write(18,*)'no em mod for grau'
c           stop 
          ENDIF
           Nlegen=max(Nlegen,NLEGENgr)
        call LEGEndre2PHASEFUNCTION(LEGENgr,NLEGENgr,
     $ 2,200,P11,ANG) 
            backgr=kextgr*salbgr*P11(2)     
c           write(*,*)'esco gr2',kextgr,salbgr,legengr(1),legengr(2),
c     $legengr(Nlegengr),legen(Nlegengr+1),Nlegen,mindex
              else
c  666  continue
           kextgr=0.0d0
           salbgr=0.0d0 
           backgr=0.0d0 
            endif
           
CCCCCCCCCCCCCC   END OF SINGLE SCATTERING PROPERTY COMPUTATIONS  CCCCCCCCCCCCCCCCCCCCCCCC   



c
ccccccc           Summing up the scattering parameters and writing the
ccccccc           input file of the scattering properties of each layer for RT3 input
c
            
            kexttot(nx,ny,nz) = kextcw +kextrr+  
     $          kextci + kextsn +  kextgr 
            kextcloud(nx,ny,nz) =max(0.0d0,kextcw)
            kextrain(nx,ny,nz) =max(0.0d0,kextrr) 
            kextice(nx,ny,nz) =max(0.0d0,kextci)
            kextsnow(nx,ny,nz) =max(0.0d0,kextsn) 
            kextgraupel(nx,ny,nz) =max(0.0d0,kextgr)
            back(nx,ny,nz)= backcw + backrr + 
     $          backci + backsn +  backgr 
c          write(18,*)'scat prop',nz,back(nx,ny,nz),
c     $ backcw,cloud_water(nx,ny,nz),backrr,rain_water(nx,ny,nz),
c     $ backci,cloud_ice(nx,ny,nz),backsn,snow(nx,ny,nz),
c     $ backgr,graupel(nx,ny,nz),backhl           
C      
        if (kexttot(nx,ny,nz).lt. 0.) write(*,*)'something wrong'       
              if (kexttot(nx,ny,nz) .le. 0. ) then
              salbtot(nx,ny,nz) = 0.0
              else
              salbtot(nx,ny,nz) = (salbcw*kextcw+  
     $ salbrr*kextrr +salbci*kextci + salbsn*kextsn+salbgr*kextgr)
     $/kexttot(nx,ny,nz)
              endif

              if ( salbtot(nx,ny,nz) .le. 0.0 ) then
              asymtot(nx,ny,nz) = 0.0d0
              else
              asymtot(nx,ny,nz) = ( asymcw*salbcw*kextcw +
     $         asymrr*salbrr*kextrr + asymci*salbci*kextci +  
     $         asymsn*salbsn*kextsn + asymgr*salbgr*kextgr)/
     $         (salbtot(nx,ny,nz)*kexttot(nx,ny,nz))
              endif
        abs(nx,ny,nz)= (1.0-salbtot(nx,ny,nz))*kexttot(nx,ny,nz)
 
        tau(nx,ny)=tau(nx,ny)+(kexttot(nx,ny,nz)+
     $ KEXTATMO(nz))*(hgt_lev(nx,ny,nz)-hgt_lev(nx,ny,nz-1))
        tau_hydro(nx,ny)=tau_hydro(nx,ny)+
     $ kexttot(nx,ny,nz)*
     $ (hgt_lev(nx,ny,nz)-hgt_lev(nx,ny,nz-1))

C   summing up the Legendre coefficient               
        
        if ( kexttot(nx,ny,nz) .le. 0.0 .or. !
     $       salbtot(nx,ny,nz) .le. 0.0) then
           FILE_PH(nz)=''
c     writing no file            

        else                    !there are hydrometeor present : a PH file is needed    

           FILE_PH(nz)='../output/tmp/PHlev'//Nzstr//'f'//frq_str ! PSG: add 'temp/'

c     write(18,*)'apri file',Nlegen,NLEGENhl,
c     $ NLEGENcw,NLEGENci,NLEGENrr,NLEGENsn,NLEGENgr
           BUIO = 31 + OMP_GET_THREAD_NUM() ! PSG: included for OpenMP
           OPEN(unit=BUIO,file =file_PH(nz), STATUS = 'unknown',
     $          form='FORMATTED') ! PSG: 31 -> BUIO and here on til close(BUIO)
           write(BUIO,*) kexttot(nx,ny,nz),'   EXINCTION'
           write(BUIO,*)  kexttot(nx,ny,nz)* salbtot(nx,ny,nz),
     $          '  SCATTERING'
           write(BUIO,*) salbtot(nx,ny,nz),'   SINGLE SCATTERING ALBEDO'
           write(BUIO,*) Nlegen-1,'      DEGREE OF LEGENDRE SERIES'
       
           do 1007 jj=1,Nlegen
          
              legen(jj)=(legencw(jj)*salbcw*kextcw +
     $          legenrr(jj)*salbrr*kextrr +  legenci(jj)*salbci*kextci +  
     $          legensn(jj)*salbsn*kextsn + legengr(jj)*salbgr*kextgr) /
     $             (salbtot(nx,ny,nz)*kexttot(nx,ny,nz)) 
           
        legen2(jj)=(legen2cw(jj)*salbcw*kextcw +
     $   legen2rr(jj)*salbrr*kextrr +    legen2ci(jj)*salbci*kextci +  
     $    legen2sn(jj)*salbsn*kextsn + legen2gr(jj)*salbgr*kextgr) /
     $      (salbtot(nx,ny,nz)*kexttot(nx,ny,nz)) 
         
        legen3(jj)=(legen3cw(jj)*salbcw*kextcw +
     $   legen3rr(jj)*salbrr*kextrr +    legen3ci(jj)*salbci*kextci +  
     $    legen3sn(jj)*salbsn*kextsn + legen3gr(jj)*salbgr*kextgr) /
     $      (salbtot(nx,ny,nz)*kexttot(nx,ny,nz))
 
          legen4(jj)=(legen4cw(jj)*salbcw*kextcw +
     $    legen4rr(jj)*salbrr*kextrr +  legen4ci(jj)*salbci*kextci +  
     $    legen4sn(jj)*salbsn*kextsn + legen4gr(jj)*salbgr*kextgr) /
     $      (salbtot(nx,ny,nz)*kexttot(nx,ny,nz)) 
   
       
          
          write(BUIO,1005)jj-1,legen(jj),legen2(jj),
     $         legen3(jj),legen4(jj),legen(jj),legen3(jj)
          g_coeff(nx,ny,nz)=legen(2)/3.0d0    ! PSG: 31 -> BUIO 
 1005     format(i3,6(1x,f9.7))    

 1007  enddo                    !end of cycle over Legendre coefficient  
       close(BUIO)    ! PSG: 31 -> BUIO 
      endif 
            
 1009 end do                    !end of cycle over the vertical layers
         
C     Preparation of the PROFILE file  (needed by RT3)
      AUIOF = 21 + OMP_GET_THREAD_NUM()        ! PSG: include Thread-dependent file unit
      OPEN(AUIOF,FILE=file_profile,FORM='FORMATTED',STATUS='unknown')   ! PSG: 21 -> AUIOF
      do 577  nz = nlyr,1,-1    ! PSG: N_lay_cut,1,-1 
         str1=''''
         offset1=index(FILE_PH(nz),' ') !position of the first blank space
         tmp_file1=FILE_PH(nz)
         FILE_PH2(nz)=str1//tmp_file1(1:offset1-1)//str1
         write(AUIOF,1013)  hgt_lev(nx,ny,nz),temp_lev(nx,ny,nz),
     $        KEXTATMO(nz), FILE_PH2(nz)   ! PSG: 21 -> AUIOF
C     write(*,*) nz, hgt_lev(nx,ny,nz),temp_lev(nx,ny,nz), ! PSG: dumn inclusion
C     $            KEXTATMO(nz)
 1013    format(f6.3,1x,f6.2,1x,E10.4,1x,a38) ! PSG: E9.4->E10.4
 577  continue                  !end of cycle over the vertical layers
      write(AUIOF,1012)  hgt_lev(nx,ny,0),temp_lev(nx,ny,0),
     $     KEXTATMO(1),''' '''   ! PSG: 21 -> AUIOF
 1012 format(f6.3,1x,f6.2,1x,E9.4,1x,a3)  
      close(AUIOF)        ! PSG: 21 -> AUIOF










!     ????? ! PSG: the following lines are commented out since micro-physics is include in netCDF

c$$$         file_profile2='output/MP/'//micro_str//'date'//date_str//  ! PSG: add 'output/'
c$$$     $'x'//xstr//'y'//ystr//'f'//frq_str 
c$$$
c$$$        OPEN(5,FILE=file_profile2,FORM='FORMATTED',STATUS='unknown')  
c$$$        write(5,1110)nx,ny, hgt_lev(nx,ny,0),
c$$$     $press_lev(nx,ny,0),temp_lev(nx,ny,0),relhum_lev(nx,ny,0)     
c$$$       do nz = 1, nlyr  
c$$$       write(5,1111)nx,ny,nz,hgt_lev(nx,ny,nz),
c$$$     $  cloud_water(nx,ny,nz),
c$$$     $  kextcloud(nx,ny,nz),rain_water(nx,ny,nz),
c$$$     $  kextrain(nx,ny,nz),   
c$$$     $ cloud_ice(nx,ny,nz),kextice(nx,ny,nz),
c$$$     $ snow(nx,ny,nz),kextsnow(nx,ny,nz),
c$$$     $ graupel(nx,ny,nz), kextgraupel(nx,ny,nz),
c$$$     $ kexttot(nx,ny,nz),salbtot(nx,ny,nz),g_coeff(nx,ny,nz),
c$$$     $ back(nx,ny,nz),KEXTATMO(nz),temp_lev(nx,ny,nz),
c$$$     $ rho_vap(nx,ny,nz)
c$$$          enddo
c$$$          close(5)
c$$$
c$$$c          write(5,*) tau(nx,ny)
c$$$ 1110    format(i3,1x,i3,1x,4(1x,f8.3))  ! PSG: ,4(1x,f7.3))
c$$$ 1111    format(i3,1x,i3,1x,i3,1x,f6.3,10(1x,e9.3),1x,e9.4,2(1x,f7.4), ! PSG: i3,1x,i3,1x,i2
c$$$     $ 1x,e9.3,1x,e9.4,1x,f5.1,1x,e9.3)

!     ????? ! PSG: end of micro-physics comment out block
  

C&&&&&&&&   I/O FILE NAMES for the MC&&&&&&&&&&&&&&&&&&
        
        
C         FILEOUT1='../output/TB/RT3TB'//date_str//micro_str//'x'//xstr
C     $        //'y'//ystr//'f'//frq_idx   ! PSG: added frq_str  (old version)
   
         FILEOUT1 = trim(NCDFOUT)//'='//date_str//'x'//xstr
     $        //'y'//ystr//'f'//frq_idx
       

      
        OUT_FILE= FILEOUT1

C        write(*,*) 'output file is: ', OUT_FILE

c$$$       NOUTLEVELS=2  ! PSG: comment out, moved to the begining of code 
c$$$       OUTLEVELS(1)=1
c$$$       OUTLEVELS(2)=nlyr+1 ! PSG: N_lay_cut+1
C       write(*,*) 'entra a '//FILE_profile//' com outlevels=',OUTLEVELS   ! PSG: test

       call  RT3(NSTOKES, NUMMU,AZIORDER, MU_VALUES,
     .          src_code,   FILE_profile, out_file,
     .          QUAD_TYPE,deltam,DIRECT_FLUX, DIRECT_MU, 
     .         GROUND_TEMP, GROUND_TYPE,
     .         GROUND_ALBEDO, GROUND_INDEX,
     .                    SKY_TEMP, WAVELENGTH, UNITS, OUTPOL,
     .                    NOUTLEVELS, OUTLEVELS,NUMAZIMUTHS,i_time)

   
       call MP_storencdf(NCDFOUT,i_time,ifreq,ny,nx,
     $      NLYR,hgt_lev(nx,ny,0:NLYR),
     $      temp_lev(nx,ny,0:NLYR), press_lev(nx,ny,0:NLYR),
     $      relhum_lev(nx,ny,1:NLYR), rho_vap(nx,ny,:),
     $      cloud_water(nx,ny,:),
     $      winddir_tmp(nx,ny,:NLYR,timeidx),
     $      windvel_tmp(nx,ny,:NLYR,timeidx),
     $      kextcloud(nx,ny,:), KEXTATMO,kexttot(nx,ny,:),
     $      salbtot(nx,ny,:), back(nx,ny,:), g_coeff(nx,ny,:))

            close(28)



 646     enddo                  ! end over time index  

    
 656  continue                ! end over ny_grid and nx_grid
!$OMP END DO
!$OMP END PARALLEL      
 777  enddo                   ! end over frequency index

      deallocate(temp_tmp, press_tmp, relhum_tmp, mixr_tmp)
      deallocate(cloud_water_tmp, rain_water_tmp, cloud_ice_tmp)
      deallocate(snow_tmp, graupel_tmp, windvel_tmp, winddir_tmp)
      deallocate(qidx)
      deallocate(TimeStamp)
          
      STOP
      END        
C     *********************************************************
C     END OF model_radtran_MW
C     *********************************************************


C     ========================================================
C     AUXILIARRY SUBROUTINES START HERE
C     -------------------------------------------------------
        subroutine get_atmosG0(temp_lev, press_lev, relhum, 
     + ngridx, ngridy, nlyr, avg_pressure, vapor_pressure,rho_vap)

c     Calculate average air pressure and vapor pressure in specified
c     layers, given the temperature, pressure, and relative humidity 
c     from cloud-resolving model output.

C+-------+---------+---------+---------+---------+---------+---------+-+
C+----------------------------------------------------------------------

      include    'parameters.inc'
c      real*8       hgt_lev(0:mxlyr)
      real*8       press_lev(mxgridx,mxgridy,0:mxlyr)
      real*8       temp_lev(mxgridx,mxgridy,0:mxlyr)
      real*8       relhum(mxgridx,mxgridy,mxlyr)
      real*8       avg_pressure(mxgridx,mxgridy,mxlyr)
      real*8       vapor_pressure(mxgridx,mxgridy,mxlyr)
      real*8      rho_vap(mxgridx,mxgridy,mxlyr)

      real*8       tavg
      real*8       esatvap

      data a0/6.107799961e0/
      data a1/4.436518521e-1/
      data a2/1.428945805e-2/
      data a3/2.650648471e-4/
      data a4/3.031240396e-6/
      data a5/2.034080948e-8/
      data a6/6.136820929e-11/

      do nx = 1, ngridx
        do ny = 1, ngridy
          do nz = 1, nlyr
          
            tavg = 0.5*(temp_lev(nx,ny,nz-1) + temp_lev(nx,ny,nz))
            avg_pressure(nx,ny,nz) = 
     +          (press_lev(nx,ny,nz) - press_lev(nx,ny,nz-1) )/
     +           log(press_lev(nx,ny,nz)/press_lev(nx,ny,nz-1))   !mb
            tc = tavg - 273.15
            ES = a0+tc*(a1+tc*(a2+tc*(a3+tc*(a4+tc*(a5+a6*tc)))))
            if (ES .lt. 0.) ES = 0.
            vapor_pressure(nx,ny,nz) =  relhum(nx,ny,nz)*ES/100.
            rho_vap(nx,ny,nz)=vapor_pressure(nx,ny,nz)*100.0
     $ /(tavg*461.5)*1e3   !g/m3
          end do
        end do
      end do  

      return
      end       




        SUBROUTINE slant_profile(Angle_view,
     $Ds_cloud_slant,Ds_obs_point,DH_bot_intersect,
     $DH_top_intersect)
C      given an observation angle, a thickness of the cloud in the
c      viewing direction (but in the horizontal) and 
c     a vector of altitudes compute the new vector of altitudes
c      INPUT=
c       Ds_cloud_slant (in km)
c       Angle_view in radiants
c      Ds_obs_point in km (positive=out of the clous, 
c        negative=within the clous)
c       hght_lev=old height vectors
c       hght_lev_new=new height vectors

      implicit none
       REAL*8 Ds_cloud_slant,Ds_obs_point,
     $ angle_view,DH_bot_intersect,DH_top_intersect

      if(Ds_obs_point.gt.0d0) then
      DH_bot_intersect=Ds_obs_point*dtan(Angle_view)
      else
      DH_bot_intersect=0.0d0
      endif
      DH_top_intersect=(Ds_obs_point+Ds_cloud_slant)*
     $dtan(Angle_view)

      RETURN
      END



          INTEGER FUNCTION LENGTH(STRING) 
cReturns length of string ignoring trailing blanks 
      CHARACTER*(*) STRING 
      DO 15, I = LEN(STRING), 1, -1 
         IF(STRING(I:I) .NE. ' ') GO TO 20 
15    CONTINUE 
20    LENGTH = I 
      END 


