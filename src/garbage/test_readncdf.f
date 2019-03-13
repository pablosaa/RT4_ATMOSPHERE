C     Test source to call the NetCDF subroutine from Fortran 77
C     this is only a simulated call from RT3/4 main code

      program test_readncdf

      include 'parameters.inc'
      character*300 filename
      character*300 pathname

      integer ngridx, ngridy, nlyr, ntime
      real*8 hgt_tmp(mxgridx,mxgridy,0:mxlyr,mxtime)
      real*8 press_tmp(mxgridx,mxgridy,0:mxlyr,mxtime)
      real*8 temp_tmp(mxgridx,mxgridy,0:mxlyr,mxtime)
      real*8 relhum_tmp(mxgridx,mxgridy,0:mxlyr,mxtime)  
      real*8 cloud_water_tmp(mxgridx,mxgridy,mxlyr,mxtime)
      real*8 rain_water_tmp(mxgridx,mxgridy,mxlyr,mxtime)
      real*8 cloud_ice_tmp(mxgridx,mxgridy,mxlyr,mxtime)
      real*8 snow_tmp(mxgridx,mxgridy,mxlyr,mxtime)
      real*8 graupel_tmp(mxgridx,mxgridy,mxlyr,mxtime)

      real*4 yy(mxtime),mm(mxtime),dd(mxtime),hh(mxtime)
      real*4 lat(mxgridx,mxgridy),lon(mxgridx,mxgridy)

      real*8 temp_lev(mxgridx,mxgridy,0:mxlyr)
      
      write(pathname,'(A)') '/home/pga082/GFI/data/RASOBS/polargmo/'
      write(filename,'(A)') 'RS_Y2013-2013_M09-09_D01-31_H00-18.nc'

      print*, 'size: ',shape(temp_tmp), ' and time:', mxtime
      call read_wyosonde(trim(pathname)//trim(filename),mxgridx,mxgridy,
     $     mxlyr,mxtime,hgt_tmp,press_tmp,temp_tmp,
     $     relhum_tmp,cloud_water_tmp,
     $     rain_water_tmp,cloud_ice_tmp,snow_tmp,graupel_tmp,
     $     ngridx,ngridy,nlyr,ntime,lat,lon,
     $     yy,mm,dd,hh)
      print*, 'saida da subroutina!'
      print*, 'ngridx', ngridx
      print*, 'ngridy', ngridy
      print*, 'nlyr', nlyr
      print*, 'ntime', ntime
      temp_lev = temp_tmp(:,:,:,ntime)
      print*,'T=[',temp_lev(1,1,0:20),'];'
      print*,'days=[',dd,'];'
      stop
      end
