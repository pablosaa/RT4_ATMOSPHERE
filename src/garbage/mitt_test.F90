module nctoys
  public
  type ncname
     character(len=15), allocatable, dimension(:) :: vars
     character(len=15), allocatable, dimension(:) :: dims
  end type ncname

  !interface getUnixTime
  !   module procedure getF2UnixTime, getCH2UnixTime
  !end interface getUnixTime
  
contains
  ! Subroutine to work with converting Date to unix-time (uses external functions)
  subroutine getF2UnixTime(datum, unixtime)
    use, intrinsic :: iso_c_binding

    implicit none
    ! Interface to the C code for Unix time retrieval:
    interface
       
       function F2UnixTime(ntime, datum, val) result(rel) bind(c, name='F2UnixTime')
         import 
         integer(kind=c_int), value :: ntime
         integer(kind=c_int) :: datum(ntime, 6)
         real(kind=c_double) :: val(ntime)
	integer(kind=c_int) :: rel
       end function F2UnixTime

    end interface
    
    integer(c_int), intent(in) :: datum(:,:)
    real(c_double), intent(out) :: unixtime(:)

    integer :: ntime
    ntime = size(datum,1)
    print*,'running on int date', ntime
	print*, datum
	!allocate(unixtime(ntime) )
	ntime = F2UnixTime(ntime, datum, unixtime)
	print*, 'passa!, ', unixtime
    return
    
  end subroutine GetF2UnixTime

  subroutine getCH2UnixTime(datum, unixtime)
    use, intrinsic :: iso_c_binding
	
    implicit none
    ! Interface to the C code for Unix time retrieval:
    interface
       
       function CH2UnixTime(ntime, nlen,  datum, val) result(res) bind(c, name='CH2UnixTime')
         import 
         integer(kind=c_int), value :: ntime, nlen
         character(kind=c_char,len=1) :: datum(ntime)
         real(kind=c_double) :: val(ntime)
	integer(kind=c_int) :: res
       end function CH2UnixTime
       
    end interface
    character(kind=c_char,len=*), intent(in) :: datum(:)
    real(kind=8), dimension(:), intent(out) :: unixtime
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
    ntime = CH2UnixTime(ntime, nlen, datum, unixtime)
	call GetF2UnixTime(date, unixtime)
	print*, 'unix is: ', unixtime
    return
    
  end subroutine GetCH2UnixTime   
 
end module


program sample
	
	use nctoys
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
	
	allocate(date(2,6), unixtt(2) )
  date(1,:) = (/2014,10,1,0,0,0/)
  date(2,:) = (/2004,5,9,0,0,0/)
!unixtt =  mitt_unixt(date)
	call getF2UnixTime(date,unixtt)

  print*, date
  print*, unixtt
	allocate(character(len=19):: datum(2) )
	datum = ["2014-10-01 00:00:00", "2004-05-09 00:00:00"]

	call getCH2UnixTime(datum, unixtt)
  stop
end program sample

