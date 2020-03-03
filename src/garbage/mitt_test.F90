module nctoys
  public
  type ncname
     character(len=15), allocatable, dimension(:) :: vars
     character(len=15), allocatable, dimension(:) :: dims
  end type ncname

  interface getUnixTime
     module procedure getF2UnixTime, getCH2UnixTime, getUnixTime2Date
  end interface getUnixTime
  !interface
  !	subroutine getUnixTime2Date(unixt, datum) bind(c, name='UnixTime2Date')
!		use, intrinsic :: iso_c_binding	
		
!		real(kind=c_double), value :: unixt
!		integer(kind=c_int) :: datum(6)
!	end subroutine getUnixTime2Date
 ! end interface
contains
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
	end function

  ! Subroutine to work with converting Date to unix-time (uses external functions)
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
    print*,'running on int date', ntime
    print*, datum
    allocate(unixtime(ntime) )
    call F2UnixTime(ntime, datum, unixtime)
    print*, 'passa!, ', ntime, unixtime
    return
    
  end function getF2UnixTime

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
    character(kind=c_char,len=*), intent(in) :: datum(:)
    real(c_double), allocatable, dimension(:) :: unixtime
    integer, allocatable, dimension(:,:) :: date
    integer :: ntime, nlen,i
    nlen = len(datum)
    ntime = size(datum,1)
    print*,'char is ', ntime, nlen
    write(*,'(A)') (datum(i), i=1,ntime)
    allocate(date(ntime,6))
    do i=1, ntime
       read(datum(i), '(I04XI02XI02XI02XI02XI02)') date(i,:)
       print*, date(i,:)
    end do
    allocate(unixtime(ntime))
    call CH2UnixTime(ntime, nlen, datum, unixtime)
    ! call GetF2UnixTime(date, unixtime)
    print*, 'unix is: ', unixtime
    return
    
  end function GetCH2UnixTime
 
end module


program sample
	
  use nctoys, only : getUnixTime
  !use, intrinsic :: iso_c_binding

  implicit none

  !interface
  !   function mitt_unixt(datum) result(val) bind(c, name='unixtime')
  !     use, intrinsic :: iso_c_binding
  !     integer(kind=c_int) :: datum(6)
  !     real(kind=c_double) :: val
  !   end function mitt_unixt
  !end interface

  !integer(c_int) :: date(6)
  !real(c_double) :: unixtt

  integer, allocatable :: date(:,:)
  real(kind=8), allocatable :: unixtt(:)
  character(len=:), allocatable :: datum(:)
	integer :: dd(6)
  character(len=8) :: datestr

  allocate(date(2,6), unixtt(2) )
  date(1,:) = (/2014,10,1,0,0,0/)
  date(2,:) = (/2004,5,9,0,0,0/)
  unixtt =  getUnixTime(date)
  !call getUnixTime(date,unixtt)

  print*, date
  print*, unixtt
  allocate(character(len=19):: datum(2) )
  datum = ["2014-10-01_00:00:00", "2004-05-09_00:00:00"]

  !call getUnixTime(datum, unixtt)
  unixtt = getUnixTime(datum)
	
	!call getUnixTime(unixtt(1), dd) !2Date(unixtt(1), dd)
	dd = getUnixTime(unixtt(1))
	print*, 'Date for unix time: ', unixtt(1), ' is ', dd
	!call getUnixTime(unixtt(2), dd) !2Date(unixtt(2), dd)
	dd = getUnixTime(unixtt(2))
	print*, 'Date for unix time: ', unixtt(2), ' is ', dd
	dd(1) = dd(1) - 2000
	write(datestr, '(4I2.2)') dd(1:4) 
	print*, 'The string is: ', datestr
  stop
end program sample

